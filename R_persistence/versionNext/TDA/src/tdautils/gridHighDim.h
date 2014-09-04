inline std::vector< unsigned char > isInternal( unsigned int argIdx, const std::vector< unsigned int > & argGridNum )
{
    std::vector< unsigned char > resIsInt;
    resIsInt.reserve( argGridNum.size() );
    std::vector< unsigned int >::const_iterator itrDim;
    for (itrDim = argGridNum.begin(); itrDim != argGridNum.end(); itrDim++)
    {
        resIsInt.push_back( static_cast<unsigned char>( argIdx % (*itrDim) > 0 ) );
        argIdx /= (*itrDim);
    }
    return resIsInt;
}

inline std::vector< std::vector< unsigned char > > verticesLessVertex( const std::vector< unsigned char > & argVtx, const bool argAlsoEqual )
{
	std::vector< std::vector< unsigned char > > resCubeVertices;
    unsigned int idxVtx, vtxNum;
    std::vector< unsigned int > oneTwoVec;
	oneTwoVec.reserve( argVtx.size() );
	std::vector< unsigned char >::const_iterator itrVtx;
	for (itrVtx = argVtx.begin(); itrVtx != argVtx.end(); ++itrVtx)
	{
		oneTwoVec.push_back( 1+static_cast<unsigned int>(*itrVtx) );
	}
	
    vtxNum = std::accumulate( oneTwoVec.begin(), oneTwoVec.end(), 1, std::multiplies< unsigned int >() );
	if ( !argAlsoEqual )
	{
		vtxNum -= 1;
	}
    resCubeVertices.reserve( vtxNum );
    for ( idxVtx = 0; idxVtx < vtxNum; ++idxVtx )
    {
        resCubeVertices.push_back( isInternal( idxVtx, oneTwoVec ) );
    }
	return resCubeVertices;
}

std::vector< std::map< std::vector< unsigned char >, std::vector< std::vector< std::vector< unsigned char > > > > > triangulateHypercube(const int argDim )
{
    std::vector< std::map< std::vector< unsigned char >, std::vector< std::vector< std::vector< unsigned char > > > > > resTriedCube;
    resTriedCube.reserve( argDim+1 );

    // vertices of hypercube
	std::vector< unsigned char > rootVtx( argDim, 1 );
    std::vector< std::vector< unsigned char > > cubeVertices = verticesLessVertex( rootVtx, true );

	std::map< std::vector< unsigned char >, std::vector< std::vector< std::vector< unsigned char > > > > mapDirSmpxVec;
	std::vector< std::vector< std::vector< unsigned char > > > dirSmpxVec;
	std::vector< std::vector< unsigned char > > dirSmpx;
	std::vector< std::vector< unsigned char > >::const_iterator itrVtx;

	// 0 dim
	mapDirSmpxVec.clear();
	for (itrVtx = cubeVertices.begin(); itrVtx != cubeVertices.end(); ++itrVtx)
	{
		dirSmpxVec.clear();
		dirSmpx.clear();
		dirSmpx.push_back( *itrVtx );
		dirSmpxVec.push_back( dirSmpx );
		mapDirSmpxVec[ *itrVtx ] = dirSmpxVec;
	}
	resTriedCube.push_back( mapDirSmpxVec );

	unsigned int idxDim;
	std::vector< std::vector< unsigned char > > vtxLessVtx;
	std::vector< std::vector< unsigned char > >::const_iterator itrLessVtx;
	std::vector< std::vector< std::vector< unsigned char > > > dirSmpxVecPrev;
	std::vector< std::vector< std::vector< unsigned char > > >::iterator itrSmpxVec;

	for (idxDim = 1; idxDim <= argDim; ++idxDim)
	{
		mapDirSmpxVec.clear();
		for (itrVtx = cubeVertices.begin(); itrVtx != cubeVertices.end(); ++itrVtx)
		{
			dirSmpxVec.clear();
			vtxLessVtx = verticesLessVertex( *itrVtx, false );
			for (itrLessVtx = vtxLessVtx.begin(); itrLessVtx != vtxLessVtx.end(); ++itrLessVtx)
			{
				dirSmpxVecPrev = resTriedCube.at( idxDim-1 ).at( *itrLessVtx );
				for ( itrSmpxVec = dirSmpxVecPrev.begin(); itrSmpxVec != dirSmpxVecPrev.end(); ++itrSmpxVec )
				{
					itrSmpxVec->push_back( *itrVtx ); 
				}
				dirSmpxVec.insert( dirSmpxVec.end(), dirSmpxVecPrev.begin(), dirSmpxVecPrev.end() );
			}
			mapDirSmpxVec[ *itrVtx ] = dirSmpxVec;
		}
		resTriedCube.push_back( mapDirSmpxVec );
	}

    return resTriedCube;
}



// add a single simplex to the filtration
inline void addSimplex(Fltr & argFltr, const double * const extFcnVal, VertexVector & argVtx)
{
    VertexVector::const_iterator itrVtx = argVtx.begin();
    double maxFcnVal = extFcnVal[ *itrVtx ];
     for (; itrVtx != argVtx.end(); ++itrVtx)
     {
        maxFcnVal = std::max( maxFcnVal, extFcnVal[ *itrVtx ] );
     }
     argFltr.push_back( Smplx(argVtx.begin(), argVtx.end(), maxFcnVal) ); 
}


void addSimplices(Fltr & argFltr, const double * const extFcnVal, const int argIdxCur, const std::vector< unsigned int > & argGridNum, const unsigned int argIdxDim, std::vector< std::map< std::vector< unsigned char >, std::vector< std::vector< std::vector< unsigned char > > > > > & argTriedCube )
{
    std::vector< unsigned char > isInt = isInternal( argIdxCur, argGridNum );
    std::vector< std::vector< std::vector< unsigned char > > > dirSmpxVec = (argTriedCube.at( argIdxDim )).at( isInt );
    std::vector< std::vector< std::vector< unsigned char > > >::const_iterator itrDirSmpxVec;
    std::vector< std::vector< unsigned char > >::const_iterator itrDirVtxVec;
	std::vector< unsigned char > diffVtx(argGridNum.size());
	std::vector< unsigned int > gridAccNum(argGridNum.size(),1);
	std::partial_sum( argGridNum.begin(), argGridNum.end()-1, gridAccNum.begin()+1, std::multiplies< unsigned int >() );

    VertexVector vtxVec;
    VertexVector::iterator itrVtxVec;
    vtxVec.resize( argIdxDim+1 );
    for (itrDirSmpxVec = dirSmpxVec.begin(); itrDirSmpxVec != dirSmpxVec.end(); ++itrDirSmpxVec)
    {
        for (itrDirVtxVec = itrDirSmpxVec->begin(), itrVtxVec = vtxVec.begin();
             itrDirVtxVec != itrDirSmpxVec->end(); 
             ++itrDirVtxVec, ++itrVtxVec)
        {
			std::transform( isInt.begin(), isInt.end(), itrDirVtxVec->begin(), diffVtx.begin(), std::minus< char >() );

            (*itrVtxVec) = ( argIdxCur - std::inner_product( gridAccNum.begin(), gridAccNum.end(), diffVtx.begin(), 0 ) );
        }
        addSimplex( argFltr, extFcnVal, vtxVec );
    }

}


int simplicesFromGridBarycenter(Fltr & argFltr, const double * const extFcnVal, const std::vector< unsigned int > & argGridNum)
{

    unsigned int idxCur, idxDim;
    const unsigned int gridNumProd = std::accumulate( argGridNum.begin(), argGridNum.end(), 1, std::multiplies< unsigned int >() );

   std::vector< std::map< std::vector< unsigned char >, std::vector< std::vector< std::vector< unsigned char > > > > >  triedCube = triangulateHypercube( argGridNum.size() );
  
//	std::vector< unsigned char > rootVtx( argGridNum.size(), 1 );
	

	

    for (idxCur = 0; idxCur < gridNumProd ; ++idxCur)
    {
        for (idxDim = 0; idxDim <= argGridNum.size(); ++idxDim)
        {
		    addSimplices(argFltr, extFcnVal, idxCur, argGridNum, idxDim, triedCube );
        }
    }
  
  return 0;
} // end simplicesFromGrid function
