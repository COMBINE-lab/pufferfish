#ifndef _ASSEMBLY_EDGE_CONSTRUCTOR_H_
#define _ASSEMBLY_EDGE_CONSTRUCTOR_H_

#include "vertexenumerator.h"

namespace TwoPaCo
{
	class AssemblyEdgeConstructor
	{
	public:
		AssemblyEdgeConstructor(const std::vector<std::string> & inputFileName, const std::string & marksFileName, const VertexEnumerator & vertexEnumerator) : 				
				vertexEnumerator_(vertexEnumerator)
		{
			int64_t vertexLength = vertexEnumerator_.GetHashSeed().VertexLength();
			int64_t edgeLength = vertexLength + 1;

			std::vector<std::vector<bool> > junctionMark;
			size_t chrNumber = 0;
			ChrReader chrReader(inputFileName);
			JunctionPositionReader junctionReader(marksFileName);
//			std::unique_ptr<ConcurrentBitVector> bloomFilter = vertexEnumerator_.ReloadBloomFilter();
			for (std::string chr; chrReader.NextChr(chr); chrNumber++)
			{								
				//Read the current vector of junction marks. Notice: should be done in order!
				junctionMark.push_back(std::vector<bool>(chr.size(), false));
				junctionReader.RestoreVector(junctionMark.back(), chrNumber);
				//Init hash function				
				VertexRollingHash hash(vertexEnumerator.GetHashSeed(), chr.begin(), vertexEnumerator.GetHashSeed().HashFunctionsNumber());
				for (int64_t i = 0; i <= int64_t(chr.size()) - edgeLength; i++)
				{
					std::string vertex = chr.substr(i, vertexLength);
					//Check if the Bloom filter contains an edge
					assert(IsOutgoingEdgeInBloomFilter(hash, *bloomFilter, chr[i + edgeLength - 1]));
					if (i > 0)
					{
						assert(IsIngoingEdgeInBloomFilter(hash, *bloomFilter, chr[i - 1]));
					}
					
					//Check the if the vertex is a junction
					if (vertexEnumerator_.GetId(vertex) != INVALID_VERTEX)
					{
						//Found a junction, check that the mark in the vector is set
						assert(junctionMark[chrNumber][i] == true);
					}

					hash.Update(chr[i], chr[i + vertexLength]);
					//Check that the hash values were updated correctly
					assert(hash.Assert(chr.begin() + i + 1));
				}
			}
		}
		
	private:
		const VertexEnumerator & vertexEnumerator_;
		DISALLOW_COPY_AND_ASSIGN(AssemblyEdgeConstructor);
	};


}

#endif