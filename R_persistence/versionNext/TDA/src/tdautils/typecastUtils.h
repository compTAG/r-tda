template<typename PersistenceDiagram, typename RcppMatrix>
inline PersistenceDiagram RcppToDionysus(const RcppMatrix& rcppMatrix) {
	PersistenceDiagram dionysusDiagram;
	const unsigned rowNum = rcppMatrix.nrow();
	for (unsigned rowIdx = 0; rowIdx < rowNum; ++rowIdx)
	{
		dionysusDiagram.push_back(typename PersistenceDiagram::Point(
				rcppMatrix[rowIdx + 0 * rowNum], rcppMatrix[rowIdx + 1 * rowNum]));
	}
	return dionysusDiagram;
}



template<typename RcppMatrix, typename StlMatrix>
inline RcppMatrix concatStlToRcpp(const std::vector< StlMatrix >& stlMatrices,
		bool includeIndex, unsigned colNum) {
	unsigned rowNum = 0;

	typename std::vector< StlMatrix >::const_iterator vecItr;
	for (vecItr = stlMatrices.begin(); vecItr != stlMatrices.end(); ++vecItr) {
		rowNum += vecItr->size();
	}
	RcppMatrix rcppMatrix(rowNum, colNum);

	unsigned vecIdx, rowIdx, colIdx;
	for (vecIdx = 0, rowIdx = 0; vecIdx < stlMatrices.size(); ++vecIdx) {
		typename StlMatrix::const_iterator matItr;
		for (matItr = stlMatrices[vecIdx].begin();
				matItr != stlMatrices[vecIdx].end(); ++matItr, ++rowIdx) {
			if (includeIndex) {
				rcppMatrix[rowIdx] = vecIdx;
				for (colIdx = 0; colIdx < colNum - 1; ++colIdx) {
					rcppMatrix[rowIdx + (colIdx + 1) * rowNum] = (*matItr)[colIdx];
				}
			}
			else {
				for (colIdx = 0; colIdx < colNum; ++colIdx) {
					rcppMatrix[rowIdx + colIdx * rowNum] = (*matItr)[colIdx];
				}
			}
		}
	}

	return rcppMatrix;
}



template<typename RcppList, typename RcppVector, typename StlSet>
inline RcppList StlToRcppList(
		const std::vector< std::vector< StlSet > >& stlSets) {
	unsigned rowNum = 0;

	typename std::vector< std::vector< StlSet > >::const_iterator vecsItr;
	for (vecsItr = stlSets.begin(); vecsItr != stlSets.end(); ++vecsItr) {
		rowNum += vecsItr->size();
	}
	RcppList rcppList(rowNum);

	typename RcppList::iterator listItr;
	typename std::vector< StlSet >::const_iterator setVecItr;
	typename StlSet::const_iterator setItr;
	unsigned setIdx;
	for (vecsItr = stlSets.begin(), listItr = rcppList.begin();
			vecsItr != stlSets.end(); ++vecsItr) {

		for (setVecItr = vecsItr->begin(); setVecItr != vecsItr->end(); ++setVecItr, ++listItr) {
			RcppVector rcppVec(setVecItr->size());
			for (setIdx = 0, setItr = setVecItr->begin(); setItr != setVecItr->end(); ++setItr, ++setIdx) {
				rcppVec[setIdx] = *setItr;
			}
			*listItr = rcppVec;
		}
	}

	return rcppList;
}