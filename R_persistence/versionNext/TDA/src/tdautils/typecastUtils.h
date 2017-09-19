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



template< typename StlMatrix, typename RcppMatrix >
inline StlMatrix RcppToStl(const RcppMatrix& rcppMatrix,
		bool is_row_names = false) {

	const unsigned rowNum = rcppMatrix.nrow();
	const unsigned colNum = rcppMatrix.ncol();
	if (is_row_names) {
		StlMatrix stlMatrix(rowNum, typename StlMatrix::value_type(colNum + 1));
		for (unsigned rowIdx = 0; rowIdx < rowNum; ++rowIdx) {
			stlMatrix[rowIdx][0] = rowIdx + 1;
		}
		for (unsigned rowIdx = 0; rowIdx < rowNum; ++rowIdx) {
			for (unsigned colIdx = 0; colIdx < colNum; ++colIdx) {
				stlMatrix[rowIdx][colIdx + 1] = rcppMatrix[rowIdx + colIdx * rowNum];
			}
		}
		return stlMatrix;
	}
	else {
		StlMatrix stlMatrix(rowNum, typename StlMatrix::value_type(colNum));
		for (unsigned rowIdx = 0; rowIdx < rowNum; ++rowIdx) {
			for (unsigned colIdx = 0; colIdx < colNum; ++colIdx) {
				stlMatrix[rowIdx][colIdx] = rcppMatrix[rowIdx + colIdx * rowNum];
			}
		}
		return stlMatrix;
	}
	
}



template< typename CGALPoint3List, typename RcppMatrix >
inline CGALPoint3List RcppToCGALPoint3(const RcppMatrix& rcppMatrix) {

	const unsigned rowNum = rcppMatrix.nrow();
  CGALPoint3List cGALPoint3List;
	for (unsigned rowIdx = 0; rowIdx < rowNum; ++rowIdx) {
    cGALPoint3List.push_back(
				typename CGALPoint3List::value_type(rcppMatrix[rowIdx],
					rcppMatrix[rowIdx + rowNum], rcppMatrix[rowIdx + 2 * rowNum]));
	}
	return cGALPoint3List;
}



template< typename CGALPointDList, typename RcppMatrix >
inline CGALPointDList RcppToCGALPointD(const RcppMatrix& rcppMatrix) {

  const unsigned rowNum = rcppMatrix.nrow();
  const unsigned colNum = rcppMatrix.ncol();
  CGALPointDList cGALPointDList;
  std::vector< double > pointD(colNum);

  for (unsigned rowIdx = 0; rowIdx < rowNum; ++rowIdx) {
    for (unsigned colIdx = 0; colIdx < colNum; ++colIdx) {
      pointD[colIdx] = rcppMatrix[rowIdx + colIdx * rowNum];
    }
    cGALPointDList.push_back(
        typename CGALPointDList::value_type(pointD.size(), pointD.begin(),
          pointD.end()));
  }
  return cGALPointDList;
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



template< typename RcppList, typename RcppMatrix, typename StlVector >
inline RcppList StlToRcppMatrixList(
	const std::vector< std::vector< std::vector< StlVector > > >& stlArrays) {
	unsigned listNum = 0;

	typename std::vector< std::vector< std::vector< StlVector > > >::const_iterator vecsItr;
	for (vecsItr = stlArrays.begin(); vecsItr != stlArrays.end(); ++vecsItr) {
		listNum += vecsItr->size();
	}
	RcppList rcppList(listNum);

	typename RcppList::iterator listItr;
	typename std::vector< std::vector< StlVector > >::const_iterator matrixItr;
	typename std::vector< StlVector >::const_iterator rowItr;
	typename StlVector::const_iterator colItr;
	unsigned rowIdx, colIdx, rowNum;
	for (vecsItr = stlArrays.begin(), listItr = rcppList.begin();
	vecsItr != stlArrays.end(); ++vecsItr) {

		for (matrixItr = vecsItr->begin(); matrixItr != vecsItr->end();
		++matrixItr, ++listItr) {
			rowNum = matrixItr->size();
			if (rowNum != 0) {
				RcppMatrix rcppMatrix(rowNum, (*matrixItr)[0].size());
				for (rowIdx = 0, rowItr = matrixItr->begin();
				rowItr != matrixItr->end(); ++rowIdx, ++rowItr) {
					for (colIdx = 0, colItr = rowItr->begin(); colItr != rowItr->end();
					++colIdx, ++colItr) {
						rcppMatrix[rowIdx + colIdx * rowNum] = *colItr;
					}
				}
				*listItr = rcppMatrix;
			}
			else {
				RcppMatrix rcppMatrix(0, 0);
				*listItr = rcppMatrix;
			}
		}
	}

	return rcppList;
}



template< typename SimplexHandle, typename SimplexTree, typename RealVector >
void filtrationGudhiOne(
    const SimplexHandle & iSpx, SimplexTree & smplxTree, const int idxShift,
    RealVector & cmplxVec, double & value, RealVector & boundaryVec) {

  const unsigned nVtx = smplxTree.dimension(iSpx) + 1;

  cmplxVec = RealVector(nVtx);
  const typename SimplexTree::Simplex_vertex_range & vtxRange =
      smplxTree.simplex_vertex_range(iSpx);
  typename RealVector::iterator iCmplxVec = cmplxVec.begin();
  for (typename SimplexTree::Simplex_vertex_iterator iVtx = vtxRange.begin();
       iVtx != vtxRange.end(); ++iVtx, ++iCmplxVec) {
    // R is 1-base, while C++ is 0-base
    *iCmplxVec = *iVtx + idxShift;
  }

  value = SimplexTree::filtration(iSpx);

  // might need to change for cubical complex
  if (nVtx > 1) {
    boundaryVec = RealVector(nVtx);
  }
  const typename SimplexTree::Boundary_simplex_range & smplxRange =
      smplxTree.boundary_simplex_range(iSpx);
  typename RealVector::iterator iBdyVec = boundaryVec.begin();
  for (typename SimplexTree::Boundary_simplex_iterator iBdySpx = 
       smplxRange.begin(); iBdySpx != smplxRange.end(); ++iBdySpx, ++iBdyVec) {
    // R is 1-base, while C++ is 0-base
    *iBdyVec = SimplexTree::key(*iBdySpx) + idxShift;
  }
}



// TODO : see whether 'const SimplexTree &' is possible
template< typename RcppList, typename RcppVector, typename SimplexTree >
inline RcppList filtrationGudhiToRcpp(SimplexTree & smplxTree) {

  const unsigned nFltr = smplxTree.num_simplices();

  RcppList cmplx(nFltr);
  RcppVector values(nFltr);
  RcppList boundary(nFltr);
  typename RcppList::iterator iCmplx = cmplx.begin();
  typename RcppVector::iterator iValue = values.begin();
  typename RcppList::iterator iBdy = boundary.begin();

  const typename SimplexTree::Filtration_simplex_range & filtration =
      smplxTree.filtration_simplex_range();

  unsigned iFill = 0;
  for (typename SimplexTree::Filtration_simplex_iterator iFltr =
       filtration.begin(); iFltr != filtration.end();
       ++iFltr, ++iCmplx, ++iValue, ++iBdy) {

    // Below two lines are only needed for computing boundary
    smplxTree.assign_key(*iFltr, iFill);
    iFill++;

    RcppVector cmplxVec;
    RcppVector boundaryVec;
    filtrationGudhiOne(*iFltr, smplxTree, 1, cmplxVec, *iValue, boundaryVec);
    *iCmplx = cmplxVec;
    *iBdy = boundaryVec;
  }

  return RcppList::create(cmplx, values, boundary);
}



template< typename SimplexTree, typename RcppVector, typename RcppList >
inline SimplexTree filtrationRcppToGudhi(const RcppList & rcppList) {

  const RcppList rcppComplex(rcppList[0]);
  const RcppVector rcppValue(rcppList[1]);
  SimplexTree smplxTree;

  typename RcppList::const_iterator iCmplx = rcppComplex.begin();
  typename RcppVector::const_iterator iValue = rcppValue.begin();
  for (; iCmplx != rcppComplex.end(); ++iCmplx, ++iValue) {
    const RcppVector rcppVec(*iCmplx);
    RcppVector gudhiVec(rcppVec.size());
    typename RcppVector::const_iterator iRcpp = rcppVec.begin();
    typename RcppVector::iterator iGudhi = gudhiVec.begin();
    for (; iRcpp != rcppVec.end(); ++iRcpp, ++iGudhi) {
      // R is 1-base, while C++ is 0-base
      *iGudhi = *iRcpp - 1;
    }
    smplxTree.insert_simplex(gudhiVec, *iValue);
  }

  return smplxTree;
}



// TODO : see whether 'const SimplexTree &' is possible
template< typename Filtration, typename SimplexTree >
inline Filtration filtrationGudhiToDionysus(SimplexTree & smplxTree) {

  const typename SimplexTree::Filtration_simplex_range & fltrGudhi =
      smplxTree.filtration_simplex_range();
  Filtration fltrDionysus;
  unsigned iFill = 0;

  for (typename SimplexTree::Filtration_simplex_iterator iSpx =
       fltrGudhi.begin(); iSpx != fltrGudhi.end(); ++iSpx) {

    // Below two lines are only needed for computing boundary
    smplxTree.assign_key(*iSpx, iFill);
    iFill++;

    std::vector< double > cmplxVec;
    double value;
    std::vector< double > boundaryVec;
    filtrationGudhiOne(*iSpx, smplxTree, 0, cmplxVec, value, boundaryVec);

    fltrDionysus.push_back(typename Filtration::Simplex(
      cmplxVec.begin(), cmplxVec.end(), value));
  }

  return fltrDionysus;
}



// TODO : see whether 'const SimplexTree &' is possible
template< typename Column, typename Dimension, typename SimplexTree,
          typename VectorList, typename RealVector, typename Boundary >
inline void filtrationGudhiToPhat(
    SimplexTree & smplxTree, VectorList & cmplx, RealVector & values,
    Boundary & boundary_matrix) {

  const unsigned nFltr = smplxTree.num_simplices();

  const typename SimplexTree::Filtration_simplex_range & fltrGudhi =
      smplxTree.filtration_simplex_range();
  unsigned iFill = 0;

  cmplx = VectorList(nFltr);
  values = RealVector(nFltr);
  boundary_matrix.set_num_cols(nFltr);
  typename VectorList::iterator iCmplx = cmplx.begin();
  typename RealVector::iterator iValue = values.begin();
  
  unsigned iCol = 0;
  for (typename SimplexTree::Filtration_simplex_iterator iSpx =
       fltrGudhi.begin(); iSpx != fltrGudhi.end();
       ++iSpx, ++iCmplx, ++iValue, ++iCol) {

    // Below two lines are only needed for computing boundary
    smplxTree.assign_key(*iSpx, iFill);
    iFill++;

    Column cmplxVec;
    Column boundary_indices;
    filtrationGudhiOne(
        *iSpx, smplxTree, 0, cmplxVec, *iValue, boundary_indices);
    *iCmplx = cmplxVec;

    std::sort(boundary_indices.begin(), boundary_indices.end());
    boundary_matrix.set_col(iCol, boundary_indices);
    Dimension dim_of_column = smplxTree.dimension(*iSpx);
    boundary_matrix.set_dim(iCol, dim_of_column);
  }
}



template< typename Simplex, typename SimplexMap, typename RealVector >
inline void filtrationDionysusOne(
  const Simplex & c, const SimplexMap & simplex_map, const int idxShift,
  RealVector & cmplxVec, double & value, RealVector & boundaryVec) {

  const unsigned nVtx = c.dimension() + 1;

  cmplxVec = RealVector(nVtx);
  typename RealVector::iterator iCmplxVec = cmplxVec.begin();
  for (typename Simplex::VertexContainer::const_iterator vit =
       c.vertices().begin(); vit != c.vertices().end(); ++vit, ++iCmplxVec) {
    // R is 1-base, while C++ is 0-base
    *iCmplxVec = *vit + idxShift;
  }

  value = c.data();

  // might need to change for cubical complex
  if (nVtx > 1) {
    boundaryVec = RealVector(nVtx);
  }
  typename RealVector::iterator iBdyVec = boundaryVec.begin();
  for (typename Simplex::BoundaryIterator bit = c.boundary_begin();
       bit != c.boundary_end(); ++bit, ++iBdyVec) {
    // R is 1-base, while C++ is 0-base
    *iBdyVec = simplex_map.find(*bit)->second + idxShift;
  }
}



template< typename RcppList, typename RcppVector, typename Filtration >
inline RcppList filtrationDionysusToRcpp(const Filtration & filtration) {

  const unsigned nFltr = filtration.size();
  std::map< typename Filtration::Simplex, unsigned,
    typename Filtration::Simplex::VertexComparison > simplex_map;
  unsigned size_of_simplex_map = 0;

  RcppList cmplx(nFltr);
  RcppVector values(nFltr);
  RcppList boundary(nFltr);
  typename RcppList::iterator iCmplx = cmplx.begin();
  typename RcppVector::iterator iValue = values.begin();
  typename RcppList::iterator iBdy = boundary.begin();

  for (typename Filtration::Index it = filtration.begin();
       it != filtration.end(); ++it, ++iCmplx, ++iValue, ++iBdy) {
    const typename Filtration::Simplex & c = filtration.simplex(it);

    RcppVector cmplxVec;
    RcppVector boundaryVec;
    filtrationDionysusOne(c, simplex_map, 1, cmplxVec, *iValue, boundaryVec);
    *iCmplx = cmplxVec;
    *iBdy = boundaryVec;

    simplex_map.insert(typename
        std::map< typename Filtration::Simplex, unsigned >::value_type(
            c, size_of_simplex_map++));
  }

  return RcppList::create(cmplx, values, boundary);
}



template< typename Filtration, typename RcppVector, typename RcppList >
inline Filtration filtrationRcppToDionysus(const RcppList & rcppList) {

  const RcppList rcppComplex(rcppList[0]);
  const RcppVector rcppValue(rcppList[1]);
  Filtration filtration;

  typename RcppList::const_iterator iCmplx = rcppComplex.begin();
  typename RcppVector::const_iterator iValue = rcppValue.begin();
  for (; iCmplx != rcppComplex.end(); ++iCmplx, ++iValue) {
    const RcppVector rcppVec(*iCmplx);
    RcppVector dionysusVec(rcppVec.size());
    typename RcppVector::const_iterator iRcpp = rcppVec.begin();
    typename RcppVector::iterator iDionysus = dionysusVec.begin();
    for (; iRcpp != rcppVec.end(); ++iRcpp, ++iDionysus) {
      // R is 1-base, while C++ is 0-base
      *iDionysus = *iRcpp - 1;
    }
    filtration.push_back(typename Filtration::Simplex(
        dionysusVec.begin(), dionysusVec.end(), *iValue));
  }

  return filtration;
}



template< typename SimplexTree, typename Filtration >
inline SimplexTree filtrationDionysusToGudhi(const Filtration & filtration) {

  std::map< typename Filtration::Simplex, unsigned,
      typename Filtration::Simplex::VertexComparison > simplex_map;
  unsigned size_of_simplex_map = 0;
  SimplexTree smplxTree;

  for (typename Filtration::Index it = filtration.begin();
       it != filtration.end(); ++it) {
    const typename Filtration::Simplex & c = filtration.simplex(it);

    std::vector< double > cmplxVec;
    double value;
    std::vector< double > boundaryVec;
    filtrationDionysusOne(c, simplex_map, 0, cmplxVec, value, boundaryVec);

    smplxTree.insert_simplex(cmplxVec, value);

    simplex_map.insert(typename std::map< typename Filtration::Simplex,
        unsigned >::value_type(c, size_of_simplex_map++));
  }

  return smplxTree;
}



template< typename Column, typename Dimension, typename Filtration,
          typename VectorList, typename RealVector, typename Boundary >
inline void filtrationDionysusToPhat(
    const Filtration & filtration, VectorList & cmplx, RealVector & values,
    Boundary & boundary_matrix) {

  const unsigned nFltr = filtration.size();
  std::map< typename Filtration::Simplex, typename Column::value_type,
    typename Filtration::Simplex::VertexComparison > simplex_map;
  typename Column::value_type size_of_simplex_map = 0;

  cmplx = VectorList(nFltr);
  values = RealVector(nFltr);
  boundary_matrix.set_num_cols(nFltr);
  typename VectorList::iterator iCmplx = cmplx.begin();
  typename RealVector::iterator iValue = values.begin();

  for (typename Filtration::Index it = filtration.begin();
       it != filtration.end(); ++it, ++iCmplx, ++iValue) {
    const typename Filtration::Simplex & c = filtration.simplex(it);

    Column cmplxVec;
    Column boundary_indices;
    filtrationDionysusOne(
        c, simplex_map, 0, cmplxVec, *iValue, boundary_indices);
    *iCmplx = cmplxVec;

    std::sort(boundary_indices.begin(), boundary_indices.end());
    boundary_matrix.set_col(size_of_simplex_map, boundary_indices);
    Dimension dim_of_column = c.dimension();
    boundary_matrix.set_dim(size_of_simplex_map, dim_of_column);

    simplex_map.insert(typename std::map< typename Filtration::Simplex,
        typename Column::value_type >::value_type(c, size_of_simplex_map++));
  }
}