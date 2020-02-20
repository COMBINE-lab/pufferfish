#ifndef _VERTEX_ENUMERATOR_H_
#define _VERTEX_ENUMERATOR_H_

#define MAX_CAPACITY 20

#include <deque>
#include <cstdio>
#include <numeric>
#include <sstream>
#include <unordered_map>

#include <tbb/tbb.h>
#include <tbb/mutex.h>
#include <tbb/compat/thread>
#include <tbb/spin_rw_mutex.h>
#include <tbb/blocked_range.h>
#include <tbb/parallel_sort.h>
#include <tbb/parallel_reduce.h>
#include <tbb/task_scheduler_init.h>
#include <tbb/concurrent_unordered_set.h>

#include "junctionapi/junctionapi.h"

#include "vertexrollinghash.h"
#include "streamfastaparser.h"
#include "bifurcationstorage.h"
#include "candidateoccurence.h"

namespace TwoPaCo
{
	class VertexEnumerator
	{
	public:
		virtual size_t GetVerticesCount() const = 0;
		virtual int64_t GetId(const std::string & vertex) const = 0;
		virtual const VertexRollingHashSeed & GetHashSeed() const = 0;
		//virtual std::unique_ptr<ConcurrentBitVector> ReloadBloomFilter() const = 0;

		virtual ~VertexEnumerator()
		{

		}
	};

	std::unique_ptr<VertexEnumerator> CreateEnumerator(const std::vector<std::string> & fileName,
		size_t vertexLength,
		size_t filterSize,
		size_t hashFunctions,
		size_t rounds,
		size_t threads,
		size_t abundance,
		const std::string & tmpFileName,
		const std::string & outFileName,
		std::ostream & logStream);

	template<size_t CAPACITY>
	class VertexEnumeratorImpl : public VertexEnumerator
	{
	private:
		std::string filterDumpFile_;
		std::vector<size_t> bufferOffset_;
		VertexRollingHashSeed hashFunctionSeed_;
		static const size_t BUF_SIZE = 1 << 24;
		BifurcationStorage<CAPACITY> bifStorage_;		
		typedef CompressedString<CAPACITY> DnaString;
		typedef CandidateOccurence<CAPACITY> Occurence;

		class FilterFillerWorker;

		class OccurenceHash
		{
		public:
			uint64_t operator()(const Occurence & occurence) const
			{
				return occurence.Hash();
			}
		};

		class DnaStringHash
		{
		public:
			uint64_t operator()(const DnaString & dnaString) const
			{
				return dnaString.Hash();
			}
		};

		class OccurenceEquality
		{
		public:
			bool operator()(const Occurence & occurence1, const Occurence & occurence2) const
			{
				return occurence1.EqualBase(occurence2);
			}
		};

		typedef tbb::concurrent_unordered_set<Occurence, OccurenceHash, OccurenceEquality> OccurenceSet;

	public:

		~VertexEnumeratorImpl<CAPACITY>()
		{
			std::remove(filterDumpFile_.c_str());
		}

		int64_t GetId(const std::string & vertex) const
		{
			return bifStorage_.GetId(vertex.begin());
		}

		size_t GetVerticesCount() const
		{
			return bifStorage_.GetDistinctVerticesCount();
		}

		const VertexRollingHashSeed & GetHashSeed() const
		{
			return hashFunctionSeed_;
		}
		/*
		std::unique_ptr<ConcurrentBitVector> ReloadBloomFilter() const
		{
			uint64_t realSize = uint64_t(1) << hashFunctionSeed_.BitsNumber();
			std::unique_ptr<ConcurrentBitVector> ret(new ConcurrentBitVector(realSize));
			ret->ReadFromFile(filterDumpFile_, false);
			return ret;
		}
		*/
		VertexEnumeratorImpl(const std::vector<std::string> & fileName,
			size_t vertexLength,
			size_t filterSize,
			size_t hashFunctions,
			size_t rounds,
			size_t threads,
			size_t abundance,
			const std::string & tmpDirName,
			const std::string & outFileNamePrefix,
			std::ostream & logStream) :
			vertexSize_(vertexLength),
			hashFunctionSeed_(hashFunctions, vertexLength, filterSize),
			filterDumpFile_(tmpDirName + "/filter.bin"),
			parsingException_("")
		{
			uint64_t realSize = uint64_t(1) << filterSize;
			logStream << "Threads = " << threads << std::endl;
			logStream << "Vertex length = " << vertexLength << std::endl;
			logStream << "Hash functions = " << hashFunctions << std::endl;
			logStream << "Filter size = " << realSize << std::endl;
			logStream << "Capacity = " << CAPACITY << std::endl;
			logStream << "Files: " << std::endl;
			for (const std::string & fn : fileName)
			{
				logStream << fn << std::endl;
			}
#ifdef LOGGING
			std::ofstream logFile((tmpDirName + "/log.txt").c_str());
			if (!logFile)
			{
				throw StreamFastaParser::Exception("Can't open the log file");
			}
#else
			std::ostream & logFile = std::cerr;
#endif

			tbb::mutex errorMutex;
			std::unique_ptr<std::runtime_error> error;
			
			size_t edgeLength = vertexLength + 1;
			std::vector<TaskQueuePtr> taskQueue(threads);
			for (size_t i = 0; i < taskQueue.size(); i++)
			{
				taskQueue[i].reset(new TaskQueue());
				taskQueue[i]->set_capacity(QUEUE_CAPACITY);
			}

			const uint64_t BIN_SIZE = max(uint64_t(1), realSize / BINS_COUNT);
			std::atomic<uint32_t> * binCounter = 0;

			bool offsetFill = false;
			size_t totalOffsets = 0;
			std::vector<size_t> offset;

			if (rounds > 1)
			{
				logStream << "Splitting the input kmers set..." << std::endl;
				std::vector<std::unique_ptr<tbb::tbb_thread> > workerThread(threads);
				binCounter = new std::atomic<uint32_t>[BINS_COUNT];
				std::fill(binCounter, binCounter + BINS_COUNT, 0);
				ConcurrentBitVector bitVector(realSize);
				for (size_t i = 0; i < workerThread.size(); i++)
				{
					InitialFilterFillerWorker worker(BIN_SIZE,
						hashFunctionSeed_,
						bitVector,
						vertexLength,
						*taskQueue[i],
						binCounter);
					workerThread[i].reset(new tbb::tbb_thread(worker));
				}

				DistributeTasks(fileName, edgeLength, taskQueue, error, errorMutex, logFile, offset, offsetFill);
				for (size_t i = 0; i < workerThread.size(); i++)
				{
					workerThread[i]->join();
				}

				if (error != 0)
				{
					throw *error;
				}
			}

			double roundSize = 0;
			if (rounds > 1)
			{
				roundSize = double(std::accumulate(binCounter, binCounter + BINS_COUNT, size_t(0))) / rounds;
			}


			logStream << std::string(80, '-') << std::endl;
			uint64_t low = 0;
			uint64_t high = realSize;
			uint64_t lowBoundary = 0;
			uint64_t totalFpCount = 0;
			uint64_t verticesCount = 0;
			std::ofstream bifurcationTempWrite((tmpDirName + "/bifurcations.bin").c_str(), ios::binary);
			if (!bifurcationTempWrite)
			{
				throw StreamFastaParser::Exception("Can't create a temp file");
			}

			time_t mark;
			offsetFill = true;
			std::vector<std::ifstream*> inMaskStorage(rounds);
			for (size_t round = 0; round < rounds; round++)
			{
				std::atomic<uint64_t> marks;
				marks = 0;
				mark = time(0);

				if (rounds > 1)
				{
					uint64_t accumulated = binCounter[lowBoundary];
					for (++lowBoundary; lowBoundary < BINS_COUNT; ++lowBoundary)
					{
						if (accumulated <= roundSize || round + 1 == rounds)
						{
							accumulated += binCounter[lowBoundary];
						}
						else
						{
							break;
						}
					}

					high = lowBoundary * BIN_SIZE;
				}
				else
				{
					high = realSize;
				}

				{
					ConcurrentBitVector bitVector(realSize);
					logStream << "Round " << round << ", " << low << ":" << high << std::endl;
					logStream << "Pass\tFilling\tFiltering" << std::endl << "1\t";
					{
						std::vector<std::unique_ptr<tbb::tbb_thread> > workerThread(threads);
						for (size_t i = 0; i < workerThread.size(); i++)
						{
							FilterFillerWorker worker(low,
								high,
								std::cref(hashFunctionSeed_),
								std::ref(bitVector),
								edgeLength,
								std::ref(*taskQueue[i]));
							workerThread[i].reset(new tbb::tbb_thread(worker));
						}

						DistributeTasks(fileName, edgeLength, taskQueue, error, errorMutex, logFile, offset, offsetFill);
						for (size_t i = 0; i < workerThread.size(); i++)
						{
							workerThread[i]->join();
						}
					}

					if (offsetFill)
					{
						offsetFill = false;
					}

					//bitVector.WriteToFile(filterDumpFile_);
					tbb::mutex maskStorageMutex;
					std::ofstream maskStorage(CandidateMaskFileName(tmpDirName, round).c_str(), ios::binary);
					if (!maskStorage)
					{
						throw std::runtime_error("Can't open a temporary file");
					}

					logStream << time(0) - mark << "\t";
					mark = time(0);
					{
						std::vector<std::unique_ptr<tbb::tbb_thread> > workerThread(threads);
						for (size_t i = 0; i < workerThread.size(); i++)
						{
							CandidateCheckingWorker worker(std::make_pair(low, high),
								hashFunctionSeed_,
								bitVector,
								vertexLength,
								*taskQueue[i],
								tmpDirName,
								marks,
								round,
								maskStorage,
								maskStorageMutex,
								error,
								errorMutex);

							workerThread[i].reset(new tbb::tbb_thread(worker));
						}

						DistributeTasks(fileName, vertexLength + 1, taskQueue, error, errorMutex, logFile, offset, offsetFill);
						for (size_t i = 0; i < taskQueue.size(); i++)
						{
							workerThread[i]->join();
						}

						if (error != 0)
						{
							throw *error;
						}
					}

					logStream << time(0) - mark << "\t" << std::endl;
				}

				mark = time(0);
				tbb::spin_rw_mutex mutex;
				logStream << "2\t";
				OccurenceSet occurenceSet(1 << 20);
				{
					tbb::mutex maskStorageMutex;
					inMaskStorage[round] = new std::ifstream(CandidateMaskFileName(tmpDirName, round).c_str(), ios::binary);
					if (!(*inMaskStorage[round]))
					{
						throw std::runtime_error("Can't open a temporary file");
					}

					std::vector<std::unique_ptr<tbb::tbb_thread> > workerThread(threads);
					for (size_t i = 0; i < workerThread.size(); i++)
					{
						CandidateFinalFilteringWorker worker(hashFunctionSeed_,
							vertexLength,
							*taskQueue[i],
							occurenceSet,
							mutex,
							tmpDirName,
							round,
							*inMaskStorage[round],
							maskStorageMutex,
							error,
							errorMutex);

						workerThread[i].reset(new tbb::tbb_thread(worker));
					}

					DistributeTasks(fileName, vertexLength + 1, taskQueue, error, errorMutex, logFile, offset, offsetFill);
					for (size_t i = 0; i < taskQueue.size(); i++)
					{
						workerThread[i]->join();
					}

					if (error != 0)
					{
						throw std::runtime_error(*error);
					}

					logStream << time(0) - mark << "\t";
				}

				mark = time(0);
				size_t falsePositives = 0;
				size_t truePositives = TrueBifurcations(occurenceSet, bifurcationTempWrite, vertexSize_, abundance, falsePositives);
				logStream << time(0) - mark << std::endl;
				logStream << "True junctions count = " << truePositives << std::endl;
				logStream << "False junctions count = " << falsePositives << std::endl;
				logStream << "Hash table size = " << occurenceSet.size() << std::endl;
				logStream << "Candidate marks count = " << marks << std::endl;
				logStream << std::string(80, '-') << std::endl;
				totalFpCount += falsePositives;
				verticesCount += truePositives;
				low = high + 1;
			}

			if (rounds > 1)
			{
				delete[] binCounter;
			}

			mark = time(0);			
			std::string bifurcationTempReadName = (tmpDirName + "/bifurcations.bin");
			bifurcationTempWrite.close();
			{
				std::ifstream bifurcationTempRead(bifurcationTempReadName.c_str(), ios::binary);
				if (!bifurcationTempRead)
				{
					throw StreamFastaParser::Exception("Can't open the temp file");
				}

				bifStorage_.Init(bifurcationTempRead, verticesCount, vertexLength, threads);
			}

			std::remove(bifurcationTempReadName.c_str());
			logStream << "Reallocating bifurcations time: " << time(0) - mark << std::endl;

			mark = time(0);			
			std::atomic<uint64_t> occurence;
			tbb::mutex currentStubVertexMutex;
			std::atomic<uint64_t> currentPiece;			
			uint64_t currentStubVertexId = verticesCount + 42;
			JunctionPositionWriter posWriter(outFileNamePrefix);
			occurence = currentPiece = 0;
			{
				tbb::mutex inMaskStorageMutex;
				std::vector<std::unique_ptr<tbb::tbb_thread> > workerThread(threads);
				for (size_t i = 0; i < workerThread.size(); i++)
				{
					EdgeConstructionWorker worker(vertexLength,
						*taskQueue[i],
						bifStorage_,
						posWriter,
						currentPiece,
						occurence,
						currentStubVertexId,
						currentStubVertexMutex,
						tmpDirName,
						inMaskStorage,
						inMaskStorageMutex,
						error,
						errorMutex);

					workerThread[i].reset(new tbb::tbb_thread(worker));
				}

				DistributeTasks(fileName, vertexLength + 1, taskQueue, error, errorMutex, logFile, offset, offsetFill);
				for (size_t i = 0; i < taskQueue.size(); i++)
				{
					workerThread[i]->join();
				}
			}

			for (size_t i = 0; i < rounds; i++)
			{
				inMaskStorage[i]->close();
				delete inMaskStorage[i];
				std::remove(CandidateMaskFileName(tmpDirName, i).c_str());
			}

			if (error != 0)
			{
				throw std::runtime_error(*error);
			}

			logStream << "True marks count: " << occurence << std::endl;
			logStream << "Edges construction time: " << time(0) - mark << std::endl;
			logStream << std::string(80, '-') << std::endl;
		}

	private:

		static const size_t QUEUE_CAPACITY = 16;
		static const uint64_t BINS_COUNT = 1 << 24;		

		static bool Within(uint64_t hvalue, uint64_t low, uint64_t high)
		{
			return hvalue >= low && hvalue <= high;
		}
/*
		static std::string CandidateMaskFileNamePart(const std::string & directory, size_t sequence, size_t pos, size_t round)
		{
			std::stringstream ss;
			ss << directory << "/" << sequence << "_" << pos << "_" << round << ".tmp";
			return ss.str();
		}
*/
		static std::string CandidateMaskFileName(const std::string & directory,  size_t round)
		{
			std::stringstream ss;
			ss << directory << "/" << "candidate" << "_" << round << ".tmp";
			return ss.str();
		}

		static void ReportError(tbb::mutex & errorMutex, std::unique_ptr<std::runtime_error> & error, const std::string & msg)
		{
			errorMutex.lock();
			if (error == 0)
			{
				error.reset(new StreamFastaParser::Exception(msg));
			}

			errorMutex.unlock();
		}

		class InitialFilterFillerWorker
		{
		public:
			InitialFilterFillerWorker(uint64_t binSize,
				const VertexRollingHashSeed & hashFunction,
				ConcurrentBitVector & filter,
				size_t vertexLength,
				TaskQueue & taskQueue,
				std::atomic<uint32_t> * binCounter) : binSize(binSize), hashFunction(hashFunction), filter(filter),
				vertexLength(vertexLength), taskQueue(taskQueue), binCounter(binCounter)
			{

			}

			void operator()()
			{
				size_t edgeLength = vertexLength + 1;
				while (true)
				{
					Task task;
					if (taskQueue.try_pop(task))
					{
						if (task.start == Task::GAME_OVER)
						{
							break;
						}

						if (task.str.size() < edgeLength)
						{
							continue;
						}

						std::vector<uint64_t> hashValue;
						size_t vertexLength = edgeLength - 1;						
						VertexRollingHash hash(hashFunction, task.str.begin(), hashFunction.HashFunctionsNumber());
						for (size_t pos = 0; pos + edgeLength - 1 < task.str.size(); ++pos)
						{							
							hashValue.clear();
							bool wasSet = true;
							char prevCh = task.str[pos];
							char nextCh = task.str[pos + edgeLength - 1];
							uint64_t startVertexHash = hash.GetVertexHash();							
							GetOutgoingEdgeHash(hash, nextCh, hashValue);
							for (auto hvalue : hashValue)
							{
								if (!filter.GetBit(hvalue))
								{
									wasSet = false;
									filter.SetBitConcurrently(hvalue);
								}
							}

							hash.Update(prevCh, nextCh);
							assert(hash.Assert(task.str.begin() + pos + 1));

							uint64_t endVertexHash = hash.GetVertexHash();
							if (!wasSet)
							{
								uint64_t value[] = { startVertexHash, endVertexHash };
								for (uint64_t v : value)
								{
									uint64_t bin = v / binSize;
									if (binCounter[bin] < MAX_COUNTER)
									{
										binCounter[bin].fetch_add(1);
									}
								}
							}
						}
					}
				}
			}

		private:
			uint64_t binSize;
			const VertexRollingHashSeed & hashFunction;
			ConcurrentBitVector & filter;
			size_t vertexLength;
			TaskQueue & taskQueue;
			std::atomic<uint32_t> * binCounter;
		};


		class CandidateCheckingWorker
		{
		public:
			CandidateCheckingWorker(std::pair<uint64_t, uint64_t> bound,
				const VertexRollingHashSeed & hashFunction,
				const ConcurrentBitVector & bitVector,
				size_t vertexLength,
				TaskQueue & taskQueue,
				const std::string & tmpDirectory,
				std::atomic<uint64_t> & marksCount,
				size_t round,
				std::ofstream & maskStorage,
				tbb::mutex & maskStorageMutex,
				std::unique_ptr<std::runtime_error> & error,				
				tbb::mutex & errorMutex) : bound(bound), hashFunction(hashFunction), bitVector(bitVector), vertexLength(vertexLength), taskQueue(taskQueue),
				tmpDirectory(tmpDirectory), marksCount(marksCount), maskStorage(maskStorage), maskStorageMutex(maskStorageMutex), error(error), errorMutex(errorMutex), round(round)
			{

			}

			void operator()()
			{
				uint64_t low = bound.first;
				uint64_t high = bound.second;
				std::vector<uint64_t> temp;
				ConcurrentBitVector candidateMask(Task::TASK_SIZE);
				while (true)
				{
					Task task;
					if (taskQueue.try_pop(task))
					{
						if (task.start == Task::GAME_OVER)
						{
							break;
						}

						if (task.str.size() < vertexLength)
						{
							continue;
						}

						size_t edgeLength = vertexLength + 1;
						if (task.str.size() >= vertexLength + 2)
						{
							candidateMask.Reset(task.str.size());
							VertexRollingHash hash(hashFunction, task.str.begin() + 1, hashFunction.HashFunctionsNumber());
							size_t definiteCount = std::count_if(task.str.begin() + 1, task.str.begin() + vertexLength + 1, DnaChar::IsDefinite);
							for (size_t pos = 1;; ++pos)
							{
								char posPrev = task.str[pos - 1];
								char posExtend = task.str[pos + vertexLength];
								assert(definiteCount == std::count_if(task.str.begin() + pos, task.str.begin() + pos + vertexLength, DnaChar::IsDefinite));
								if (Within(hash.GetVertexHash(), low, high) && definiteCount == vertexLength)
								{
									size_t inCount = DnaChar::IsDefinite(posPrev) ? 0 : 2;
									size_t outCount = DnaChar::IsDefinite(posExtend) ? 0 : 2;
									for (int i = 0; i < DnaChar::LITERAL.size() && inCount < 2 && outCount < 2; i++)
									{
										char nextCh = DnaChar::LITERAL[i];
										if (nextCh == posPrev || IsIngoingEdgeInBloomFilter(hash, bitVector, nextCh))
										{
											++inCount;
										}

										if (nextCh == posExtend || IsOutgoingEdgeInBloomFilter(hash, bitVector, nextCh))
										{
											++outCount;
										}
									}

									if (inCount > 1 || outCount > 1)
									{
										++marksCount;
										candidateMask.SetBitConcurrently(pos);
									}
								}

								if (pos + edgeLength < task.str.size())
								{
									char posPrev = task.str[pos];
									definiteCount += (DnaChar::IsDefinite(task.str[pos + vertexLength]) ? 1 : 0) - (DnaChar::IsDefinite(task.str[pos]) ? 1 : 0);
									hash.Update(posPrev, posExtend);
									assert(hash.Assert(task.str.begin() + pos + 1));
								}
								else
								{
									break;
								}
							}

							try
							{
								maskStorageMutex.lock();
								candidateMask.WriteToFile(maskStorage, task.offset, task.str.size());
								maskStorageMutex.unlock();
							}
							catch (std::runtime_error & err)
							{
								ReportError(errorMutex, error, err.what());
							}
						}
					}
				}
			}

		private:
			std::pair<uint64_t, uint64_t> bound;
			const VertexRollingHashSeed & hashFunction;
			const ConcurrentBitVector & bitVector;
			size_t vertexLength;
			TaskQueue & taskQueue;
			const std::string & tmpDirectory;
			std::atomic<uint64_t> & marksCount;
			size_t round;
			std::ofstream & maskStorage;
			tbb::mutex & maskStorageMutex;
			std::unique_ptr<std::runtime_error> & error;
			tbb::mutex & errorMutex;
		};



		class CandidateFinalFilteringWorker
		{
		public:
			CandidateFinalFilteringWorker(const VertexRollingHashSeed & hashFunction,
				size_t vertexLength,
				TaskQueue & taskQueue,
				OccurenceSet & occurenceSet,
				tbb::spin_rw_mutex & mutex,
				const std::string & tmpDirectory,
				size_t round,
				std::ifstream & maskStorage,
				tbb::mutex & maskStorageMutex,
				std::unique_ptr<std::runtime_error> & error,
				tbb::mutex & errorMutex) : hashFunction(hashFunction), vertexLength(vertexLength), taskQueue(taskQueue), occurenceSet(occurenceSet),
				mutex(mutex), tmpDirectory(tmpDirectory), round(round), maskStorage(maskStorage), maskStorageMutex(maskStorageMutex), error(error), errorMutex(errorMutex)
			{

			}

			void operator()()
			{
				ConcurrentBitVector candidateMask(Task::TASK_SIZE);
				while (true)
				{
					Task task;
					if (taskQueue.try_pop(task))
					{
						if (task.start == Task::GAME_OVER)
						{
							break;
						}

						if (task.str.size() < vertexLength)
						{
							continue;
						}

						size_t edgeLength = vertexLength + 1;
						if (task.str.size() >= vertexLength + 2)
						{
							VertexRollingHash hash(hashFunction, task.str.begin() + 1, 1);
							{
								try
								{
									maskStorageMutex.lock();							
									candidateMask.ReadFromFile(maskStorage, task.offset, task.str.size());
									maskStorageMutex.unlock();
								}
								catch (std::runtime_error & err)
								{
									ReportError(errorMutex, error, err.what());
								}
							}

							for (size_t pos = 1;; ++pos)
							{
								char posPrev = task.str[pos - 1];
								char posExtend = task.str[pos + vertexLength];
								if (candidateMask.GetBit(pos))
								{
									Occurence now;
									bool isBifurcation = false;						
									now.Set(hash.RawPositiveHash(0),
										hash.RawNegativeHash(0),
										task.str.begin() + pos,
										vertexLength,
										posExtend,
										posPrev,
										isBifurcation);
									size_t inUnknownCount = now.Prev() == 'N' ? 1 : 0;
									size_t outUnknownCount = now.Next() == 'N' ? 1 : 0;
									auto ret = occurenceSet.insert(now);
									typename OccurenceSet::iterator it = ret.first;
									it->Inc();
									if (!ret.second && !it->IsBifurcation())
									{
										inUnknownCount += DnaChar::IsDefinite(it->Prev()) ? 0 : 1;
										outUnknownCount += DnaChar::IsDefinite(it->Next()) ? 0 : 1;
										if (isBifurcation || it->Next() != now.Next() || it->Prev() != now.Prev() || inUnknownCount > 1 || outUnknownCount > 1)
										{
											it->MakeBifurcation();
										}
									}
								}

								if (pos + edgeLength < task.str.size())
								{
									char posPrev = task.str[pos];
									hash.Update(posPrev, posExtend);
									assert(hash.Assert(task.str.begin() + pos + 1));
								}
								else
								{
									break;
								}
							}
						}
					}
				}
			}

		private:
			const VertexRollingHashSeed & hashFunction;
			size_t vertexLength;
			TaskQueue & taskQueue;
			OccurenceSet & occurenceSet;
			tbb::spin_rw_mutex & mutex;
			const std::string & tmpDirectory;
			size_t round;
			std::ifstream & maskStorage;
			tbb::mutex & maskStorageMutex;
			std::unique_ptr<std::runtime_error> & error;
			tbb::mutex & errorMutex;
		};

		struct EdgeResult
		{
			uint32_t pieceId;
			std::vector<JunctionPosition> junction;
		};

		static bool FlushEdgeResults(std::deque<EdgeResult> & result,
			JunctionPositionWriter & writer,
			std::atomic<uint64_t> & currentPiece)
		{
			if (result.size() > 0 && result.front().pieceId == currentPiece)
			{
				for (auto junction : result.front().junction)
				{
					writer.WriteJunction(junction);
				}

				++currentPiece;
				result.pop_front();
				return true;
			}

			return false;
		}

		class EdgeConstructionWorker
		{
		public:
			EdgeConstructionWorker(size_t vertexLength,
				TaskQueue & taskQueue,
				const BifurcationStorage<CAPACITY> & bifStorage,
				JunctionPositionWriter & writer,
				std::atomic<uint64_t> & currentPiece,
				std::atomic<uint64_t> & occurences,
				uint64_t & currentStubVertexId,
				tbb::mutex & currentStubVertexMutex,
				const std::string & tmpDirectory,
				std::vector<std::ifstream*> & inMaskStorage,
				tbb::mutex & inMaskStorageMutex,
				std::unique_ptr<std::runtime_error> & error,
				tbb::mutex & errorMutex) : vertexLength(vertexLength), taskQueue(taskQueue), bifStorage(bifStorage),
				writer(writer), currentPiece(currentPiece), occurences(occurences), tmpDirectory(tmpDirectory),
				error(error), errorMutex(errorMutex), currentStubVertexId(currentStubVertexId), currentStubVertexMutex(currentStubVertexMutex), 
				inMaskStorage(inMaskStorage), inMaskStorageMutex(inMaskStorageMutex)
			{

			}							

			void operator()()
			{
				try
				{
					DnaString bitBuf;
					std::deque<EdgeResult> result;
					ConcurrentBitVector temporaryMask(Task::TASK_SIZE);
					ConcurrentBitVector candidateMask(Task::TASK_SIZE);
					while (true)
					{
						Task task;
						if (taskQueue.try_pop(task))
						{
							if (task.start == Task::GAME_OVER)
							{
								break;
							}

							if (task.str.size() < vertexLength)
							{
								continue;
							}

							size_t edgeLength = vertexLength + 1;
							if (task.str.size() >= vertexLength + 2)
							{																
								try
								{	
									candidateMask.Reset(task.str.size());
									inMaskStorageMutex.lock();
									for (size_t i = 0; i < inMaskStorage.size(); i++)
									{
										temporaryMask.ReadFromFile(*inMaskStorage[i], task.offset, task.str.size());
										candidateMask.MergeOr(temporaryMask, task.str.size());
									}

									inMaskStorageMutex.unlock();
								}
								catch (std::runtime_error & err)
								{
									ReportError(errorMutex, error, err.what());
								}

								size_t taskstart = task.start;

								EdgeResult currentResult;
								currentResult.pieceId = task.piece;
								size_t definiteCount = std::count_if(task.str.begin() + 1, task.str.begin() + vertexLength + 1, DnaChar::IsDefinite);
								for (size_t pos = 1;; ++pos)
								{
									while (result.size() > 0 && FlushEdgeResults(result, writer, currentPiece));
									int64_t bifId(INVALID_VERTEX);
									assert(definiteCount == std::count_if(task.str.begin() + pos, task.str.begin() + pos + vertexLength, DnaChar::IsDefinite));
									if (definiteCount == vertexLength && candidateMask.GetBit(pos))
									{
										bifId = bifStorage.GetId(task.str.begin() + pos);
										if (bifId != INVALID_VERTEX)
										{
											occurences++;
											currentResult.junction.push_back(JunctionPosition(task.seqId, task.start + pos - 1, bifId));
										}
									}

									if (((task.start == 0 && pos == 1) || (task.isFinal && pos == task.str.size() - vertexLength - 1)) && bifId == INVALID_VERTEX)
									{
										occurences++;
										currentStubVertexMutex.lock();								
										currentResult.junction.push_back(JunctionPosition(task.seqId, task.start + pos - 1, currentStubVertexId++));
										currentStubVertexMutex.unlock();
									}

									if (pos + edgeLength < task.str.size())
									{
										definiteCount += (DnaChar::IsDefinite(task.str[pos + vertexLength]) ? 1 : 0) - (DnaChar::IsDefinite(task.str[pos]) ? 1 : 0);
									}
									else
									{
										break;
									}
								}

								result.push_back(currentResult);
							}
						}
					}

					while (result.size() > 0)
					{
						FlushEdgeResults(result, writer, currentPiece);
					}
				}
				catch (std::runtime_error & e)
				{
					errorMutex.lock();
					error.reset(new std::runtime_error(e));
					errorMutex.unlock();
				}
			}

		private:
			size_t abundance;
			size_t vertexLength;
			TaskQueue & taskQueue;
			uint64_t & currentStubVertexId;
			const BifurcationStorage<CAPACITY> & bifStorage;
			JunctionPositionWriter & writer;
			std::atomic<uint64_t> & currentPiece;
			std::atomic<uint64_t> & occurences;
			const std::string & tmpDirectory;
			std::unique_ptr<std::runtime_error> & error;
			tbb::mutex & inMaskStorageMutex;
			std::vector<std::ifstream*> & inMaskStorage;
			tbb::mutex & errorMutex;
			tbb::mutex & currentStubVertexMutex;
		};		
		
		class FilterFillerWorker
		{
		public:
			FilterFillerWorker(uint64_t low,
				uint64_t high,
				const VertexRollingHashSeed & hashFunction,
				ConcurrentBitVector & filter,
				size_t edgeLength,
				TaskQueue & taskQueue) : low(low), high(high), hashFunction(hashFunction), filter(filter), edgeLength(edgeLength), taskQueue(taskQueue)
			{

			}

			void operator()()
			{
				std::vector<uint64_t> setup;
				std::vector<uint64_t> hashValue;
				const char DUMMY_CHAR = DnaChar::LITERAL[0];
				const char REV_DUMMY_CHAR = DnaChar::ReverseChar(DUMMY_CHAR);
				while (true)
				{
					Task task;
					if (taskQueue.try_pop(task))
					{
						if (task.start == Task::GAME_OVER)
						{
							break;
						}

						if (task.str.size() < edgeLength)
						{
							continue;
						}


						uint64_t fistMinHash0;
						uint64_t secondMinHash0;
						size_t vertexLength = edgeLength - 1;
						size_t definiteCount = std::count_if(task.str.begin(), task.str.begin() + vertexLength, DnaChar::IsDefinite);
						VertexRollingHash hash(hashFunction, task.str.begin(), hashFunction.HashFunctionsNumber());						
						for (size_t pos = 0;; ++pos)
						{
							hashValue.clear();
							char prevCh = task.str[pos];
							char nextCh = task.str[pos + edgeLength - 1];
							assert(definiteCount == std::count_if(task.str.begin() + pos, task.str.begin() + pos + vertexLength, DnaChar::IsDefinite));
							if (definiteCount == vertexLength)
							{
								fistMinHash0 = hash.GetVertexHash();
								if (DnaChar::IsDefinite(nextCh))
								{
									GetOutgoingEdgeHash(hash, nextCh, hashValue);
								}
								else
								{
									GetOutgoingEdgeHash(hash, DUMMY_CHAR, hashValue);
									GetOutgoingEdgeHash(hash, REV_DUMMY_CHAR, hashValue);
								}

								if (pos > 0 && !DnaChar::IsDefinite(task.str[pos - 1]))
								{
									GetIngoingEdgeHash(hash, DUMMY_CHAR, hashValue);
									GetIngoingEdgeHash(hash, REV_DUMMY_CHAR, hashValue);
								}
							}

							hash.Update(prevCh, nextCh);
							assert(hash.Assert(task.str.begin() + pos + 1));
							if (definiteCount == vertexLength)
							{
								secondMinHash0 = hash.GetVertexHash();
								if (Within(fistMinHash0, low, high) || Within(secondMinHash0, low, high))
								{
									for (uint64_t value : hashValue)
									{
										setup.push_back(value);
									}
								}
							}

							if (pos + vertexLength < task.str.size() - 1)
							{
								definiteCount += (DnaChar::IsDefinite(task.str[pos + vertexLength]) ? 1 : 0) - (DnaChar::IsDefinite(prevCh) ? 1 : 0);
							}
							else
							{
								break;
							}
						}
					}

					for (uint64_t hashValue : setup)
					{
						if (!filter.GetBit(hashValue))
						{
							filter.SetBitConcurrently(hashValue);
						}
					}

					setup.clear();
				}
			}

		private:
			uint64_t low;
			uint64_t high;
			const VertexRollingHashSeed & hashFunction;
			ConcurrentBitVector & filter;
			size_t edgeLength;
			TaskQueue & taskQueue;
		};


		static void DistributeTasks(const std::vector<std::string> & fileName,
			size_t overlapSize,
			std::vector<TaskQueuePtr> & taskQueue,
			std::unique_ptr<std::runtime_error> & error,
			tbb::mutex & errorMutex,
			std::ostream & logFile,
			std::vector<size_t> & offset,
			bool offsetFill)
		{
			size_t record = 0;
			size_t nowQueue = 0;
			uint32_t pieceCount = 0;
			bool fail = false;
#ifdef LOGGING
			logFile << "Starting a new stage" << std::endl;
#endif
			if (offsetFill)
			{
				offset.push_back(0);
			}

			for (size_t file = 0; file < fileName.size(); file++)
			{
#ifdef LOGGING
				logFile << "Reading " << fileName[file] << std::endl;
#endif
				const std::string & nowFileName = fileName[file];
				for (StreamFastaParser parser(nowFileName); !fail && parser.ReadRecord(); record++)
				{
					{
						errorMutex.lock();
						if (error != 0)
						{
							throw *error;
						}

						errorMutex.unlock();
					}

					std::stringstream ss;
#ifdef LOGGING
					logFile << "Processing sequence " << parser.GetCurrentHeader() << " " << ss.str() << std::endl;
#endif
					char ch;
					uint64_t prev = 0;
					uint64_t start = 0;
					std::string buf = "N";
					bool over = false;
					do
					{
						try
						{
							over = !parser.GetChar(ch);
						}
						catch (const StreamFastaParser::Exception & e)
						{
							errorMutex.lock();
							error.reset(new std::runtime_error(e.what()));
							errorMutex.unlock();
							fail = over = true;
							break;
						}

						if (!over)
						{
							start++;
							buf.push_back(DnaChar::IsDefinite(ch) ? ch : 'N');
						}

						if (buf.size() >= overlapSize && (buf.size() == Task::TASK_SIZE || over))
						{						
							for (bool found = false; !found; nowQueue = nowQueue + 1 < taskQueue.size() ? nowQueue + 1 : 0)
							{
								TaskQueuePtr & q = taskQueue[nowQueue];
								if (q->capacity() - q->size() > 0)
								{
									std::string overlap;
									if (!over)
									{
										overlap.assign(buf.end() - overlapSize, buf.end());
									}
									else
									{
										buf.push_back('N');
									}

									size_t currentTaskSize = buf.size();
									uint64_t currentOffset = offset.size() > 0 ? offset[pieceCount] : 0;
									q->push(Task(record, prev, pieceCount++, currentOffset, over, std::move(buf)));

									if (offsetFill)
									{
										offset.push_back(offset.back() + currentTaskSize / 32 + 1);
									}

#ifdef LOGGING
									logFile << "Passed chunk " << prev << " to worker " << nowQueue << std::endl;
#endif
									prev = start - overlapSize + 1;
									buf.swap(overlap);
									found = true;
								}
							}

						}

					} while (!over);
				}
			}

			for (size_t i = 0; i < taskQueue.size(); i++)
			{
				TaskQueuePtr & q = taskQueue[nowQueue];
				while (!taskQueue[i]->try_push(Task(0, Task::GAME_OVER, 0, 0, true, std::string())))
				{
					
				}
			}
		}

		uint64_t TrueBifurcations(const OccurenceSet & occurenceSet, std::ofstream & out, size_t vertexSize, size_t abundance, size_t & falsePositives) const
		{
			uint64_t truePositives = falsePositives = 0;
			for (auto it = occurenceSet.begin(); it != occurenceSet.end();++it)
			{
				bool bifurcation = it->IsBifurcation();
				if (bifurcation && it->Count() <= abundance)
				{
					++truePositives;
					it->GetBase().WriteToFile(out);
					if (!out)
					{
						throw StreamFastaParser::Exception("Can't write to a temporary file");
					}
				}
				else
				{
					falsePositives++;
				}
			}

			return truePositives;
		}

		size_t vertexSize_;
		StreamFastaParser::Exception parsingException_;
		DISALLOW_COPY_AND_ASSIGN(VertexEnumeratorImpl<CAPACITY>);
	};	
}

#endif	
