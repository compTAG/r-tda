#include <utilities/log.h>

#include <topology/simplex.h>
#include <topology/filtration.h>
#include <topology/static-persistence.h>
#include <topology/dynamic-persistence.h>
#include <topology/persistence-diagram.h>
#include <utilities/indirect.h>

#include <vector>
#include <map>
#include <iostream>


#if 1
#include <fstream>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#endif


typedef         unsigned                                            Vertex;
typedef         Simplex<Vertex, double>                             Smplx;
typedef         Smplx::VertexContainer				    VertexCont;
typedef         std::vector<Vertex>                                 VertexVector;
typedef         Filtration<Smplx>                                   Fltr;
typedef         StaticPersistence<>                                 Persistence;
typedef         PersistenceDiagram<>                                PDgm;
typedef         OffsetBeginMap<Persistence, Fltr, 
                               Persistence::iterator, 
                               Fltr::Index>                         PersistenceFiltrationMap;
typedef         OffsetBeginMap<Fltr, Persistence,
                               Fltr::Index, 
                               Persistence::iterator>               FiltrationPersistenceMap;


// add a single edge to the filtration
void addEdge(Fltr& filtr, const std::vector<double> fcnvalues, 
            int vert01, int vert02)
{
     VertexVector vertices(3);
     vertices[0] = vert01;
     vertices[1] = vert02;
     VertexVector::const_iterator bg = vertices.begin();

     double value = std::max(fcnvalues.at(vert01), fcnvalues.at(vert02));
     filtr.push_back(Smplx(bg, bg + 2, value)); 
         // std::max(fcnvalues.at(vert03),fcnvalues.at(vert04))),
} // end function to add a single edge

// add a single triangle to the filtration
void addTri(Fltr& filtr, const std::vector<double> fcnvalues, 
            int vert01, int vert02, int vert03)
{
     VertexVector vertices(3);
     vertices[0] = vert01;
     vertices[1] = vert02;
     vertices[2] = vert03;
     VertexVector::const_iterator bg = vertices.begin();

     double value = std::max(std::max(fcnvalues.at(vert01), fcnvalues.at(vert02)),
                    fcnvalues.at(vert03));
     filtr.push_back(Smplx(bg, bg + 3, value)); 
         // std::max(fcnvalues.at(vert03),fcnvalues.at(vert04))),
} // end function to add a single triangle

// add a single tet to the filtration
void addTet(Fltr& filtr, const std::vector<double> fcnvalues, 
            int vert01, int vert02, int vert03, int vert04)
{
     VertexVector vertices(3);
     vertices[0] = vert01;
     vertices[1] = vert02;
     vertices[2] = vert03;
     vertices[3] = vert04;
     VertexVector::const_iterator bg = vertices.begin();

     double value = std::max(std::max(fcnvalues.at(vert01), fcnvalues.at(vert02)),
                    std::max(fcnvalues.at(vert03), fcnvalues.at(vert04)));
     filtr.push_back(Smplx(bg, bg + 4, value)); 
         // std::max(fcnvalues.at(vert03),fcnvalues.at(vert04))),
} // end function to add a single tet

void addAllEdges(Fltr& filtr, const std::vector<double> fcnvalues, 
              const int ncols, const int nrows, int i, int j, int k)
{     
     int curidx = i + ncols*j + ncols*nrows*k;

     // ... add edge (i-1,j,k) <--> (i,j,k)
     if (i > 0)
     {
        addEdge(filtr, fcnvalues, curidx, curidx -1);
     }
     
     // ... add edge (i,j-1,k) <--> (i,j,k)
     if (j > 0)
     {
        addEdge(filtr, fcnvalues, curidx, curidx - ncols);
     }
   
     // ... add edge (i,j,k-1) <--> (i,j,k)
     if (k > 0)
     {
        addEdge(filtr, fcnvalues, curidx, curidx - nrows*ncols);
     }


     // TODO: add the rest of the code for creating edges to here
     //       from fcn simplicesFromGrid 

     // ... consider two cases for the cubical decomposition:
     if ((i+j+k)%2 == 0)
     {
	// ... EVEN BOX 
        if (i > 0 && j > 0) // top
        { 
           addEdge(filtr, fcnvalues, curidx, curidx - ncols -1);
        }
        if (i > 0 && k > 0) // back
        {
           addEdge(filtr, fcnvalues, curidx, curidx - nrows*ncols -1);
        }
        if (j > 0 && k > 0) // right
        {
           addEdge(filtr, fcnvalues, curidx, curidx - nrows*ncols - ncols);
        }
     }
     else
     {
     	// ... ODD BOX
        if (i > 0 && j > 0) // top
        {
           addEdge(filtr, fcnvalues, curidx - 1, curidx - ncols);
        }
        if (i > 0 && k > 0) // back
        {
           addEdge(filtr, fcnvalues, curidx - 1, curidx - nrows*ncols);
        }
        if (j > 0 && k > 0) // right
        {
           addEdge(filtr, fcnvalues, curidx - ncols, curidx - nrows*ncols);
        }
     }


     return;
} // end function addEdges

void addEvenTets(Fltr& filtr, const std::vector<double> fcnvalues, 
                  const int ncols, const int nrows, int i, int j, int k)
{
     assert(i > 0 && j > 0 && k > 0);
     int curidx = i + ncols*j + ncols*nrows*k;
     
     // top vertex (i, j-1, k)
     addTet(filtr, fcnvalues, curidx, curidx - 1 - ncols, curidx - ncols - nrows*ncols, curidx - ncols);
     
     // top vertex (i-1, j, k)
     addTet(filtr, fcnvalues, curidx, curidx - 1, curidx - nrows*ncols - 1, curidx -1 - ncols);

     // top vertex (i, j, k-1)
     addTet(filtr, fcnvalues, curidx, curidx - 1 - nrows*ncols, curidx - ncols - nrows*ncols, curidx - nrows*ncols);

     // top vertex (i-1, j-1, k-1)
     addTet(filtr, fcnvalues, curidx - 1 - nrows*ncols, curidx - ncols - nrows*ncols, curidx - 1 - ncols, curidx - 1 - ncols - nrows*ncols);
     
     return;
} // end fcn to add four EVEN tets

void addOddTets(Fltr& filtr, const std::vector<double> fcnvalues, 
                  const int ncols, const int nrows, int i, int j, int k)
{
     assert(i > 0 && j > 0 && k > 0);
     int curidx = i + ncols*j + ncols*nrows*k;
      
     VertexVector vertices(4);
     vertices[0] = curidx;  vertices[3] = curidx;
     vertices[1] = -1; vertices[2] = -1; 
     VertexVector::const_iterator bg = vertices.begin();
     VertexVector::const_iterator end = vertices.end();
    
     int v1, v2, v3, v4;
     double value, value2;  // max of value and value 2 is the fcn value. 

     // top vertex (i, j, k)
     v1 = curidx -1;   vertices[0] = v1;
     v2 = curidx - ncols;  vertices[1] = v2;  
     value = std::max(fcnvalues.at(v1),fcnvalues.at(v2));
     
     v3 = curidx - nrows*ncols;  vertices[2] = v3;
     v4 = curidx; vertices[3] = v4;
     value2 = std::max(fcnvalues.at(v3),fcnvalues.at(v4));
     
     filtr.push_back(Smplx(bg, bg + 4, std::max(value,value2)));
     
     // top vertex (i-1, j-1, k)
     v3 = curidx - 1 - ncols - nrows*ncols;  vertices[2] = v3;
     v4 = curidx - 1 -ncols; vertices[3] = v4;
     value2 = std::max(fcnvalues.at(v3),fcnvalues.at(v4));
     
     filtr.push_back(Smplx(bg, bg + 4, std::max(value,value2)));

     // top vertex (i, j-1, k-1)
     v1 = curidx - nrows*ncols;   vertices[0] = v1;
     v2 = curidx - 1 - ncols - nrows*ncols;  vertices[1] = v2;  
     value = std::max(fcnvalues.at(v1),fcnvalues.at(v2));
     
     v3 = curidx - ncols;  vertices[2] = v3;
     v4 = curidx -ncols - nrows*ncols; vertices[3] = v4;
     value2 = std::max(fcnvalues.at(v3),fcnvalues.at(v4));
     
     filtr.push_back(Smplx(bg, bg + 4, std::max(value,value2)));

     // top vertex (i-1, j, k-1)
     v3 = curidx - 1;  vertices[2] = v3;
     v4 = curidx -1 - nrows * ncols; vertices[3] = v4;
     value2 = std::max(fcnvalues.at(v3),fcnvalues.at(v4));
    
     filtr.push_back(Smplx(bg, bg + 4, std::max(value,value2)));
      

     return;
} // end fcn addEvenTets

void addTriNTet(Fltr& filtr, const std::vector<double> fcnvalues, 
                  const int ncols, const int nrows, int i, int j, int k)
{
     int curidx = i + ncols*j + ncols*nrows*k;
     
     // ... consider two cases for the cubical decomposition:
     if ((i+j+k)%2 == 0)
     {
	// ... EVEN BOX
        if (i > 0 && j > 0) // top
        {
     
           addTri(filtr, fcnvalues, curidx, curidx - ncols - 1, curidx - ncols);
           addTri(filtr, fcnvalues, curidx, curidx -1, curidx - ncols -1); 
        }
        if (i > 0 && k > 0) // back
        {
           addTri(filtr, fcnvalues, curidx, curidx - nrows*ncols -1, curidx - 1);
           addTri(filtr, fcnvalues, curidx, curidx - nrows*ncols, curidx - nrows*ncols -1);
        }
        
        if (j > 0 && k > 0) // right
        {
           addTri(filtr, fcnvalues, curidx, curidx - nrows*ncols - ncols, curidx - nrows*ncols);
           addTri(filtr, fcnvalues, curidx, curidx - ncols, curidx - nrows*ncols - ncols);
           
           if (i > 0) // middle
           {
              addTri(filtr, fcnvalues, curidx, curidx - ncols -1, curidx - ncols - nrows*ncols);
              addTri(filtr, fcnvalues, curidx, curidx -1 -nrows*ncols, curidx - ncols -1);
              addTri(filtr, fcnvalues, curidx -1 - nrows*ncols, curidx - ncols - nrows*ncols, curidx);
              addTri(filtr, fcnvalues, curidx -1 - nrows*ncols, curidx - 1 - ncols, curidx - ncols - nrows*ncols);
             
              // ... add center tets 
              addTet(filtr, fcnvalues, curidx -1 - nrows*ncols, curidx - 1 - ncols, curidx - ncols - nrows*ncols, curidx);
              // ... add remaining tets 
              addEvenTets(filtr, fcnvalues, ncols, nrows, i, j, k);
           }
        }
     } // end if for even case 
     else {
        // ... ODD CASE
        if (i > 0 && j > 0) // top
        {
           addTri(filtr, fcnvalues, curidx -1, curidx - ncols, curidx);
           addTri(filtr, fcnvalues, curidx -1, curidx - ncols -1, curidx - ncols);
        }

        if (i > 0 && k > 0) // back
        {
           addTri(filtr, fcnvalues, curidx -1, curidx - nrows*ncols, curidx - nrows*ncols - 1);
           addTri(filtr, fcnvalues, curidx -1, curidx, curidx - nrows*ncols);
        }  
        
        if (j > 0 && k > 0) // right
        {
           addTri(filtr, fcnvalues, curidx - ncols, curidx - nrows*ncols, curidx - ncols - nrows*ncols);
           addTri(filtr, fcnvalues, curidx - ncols, curidx, curidx - nrows*ncols);
           
           if ( i > 0) // middle
           { 
              addTri(filtr, fcnvalues, curidx -1, curidx - ncols, curidx - nrows*ncols);
              addTri(filtr, fcnvalues, curidx -1, curidx - nrows*ncols - ncols - 1, curidx - ncols);
              addTri(filtr, fcnvalues, curidx - nrows*ncols, curidx - nrows*ncols - ncols -1, curidx - ncols);
              addTri(filtr, fcnvalues, curidx - nrows*ncols, curidx - 1, curidx - nrows*ncols - ncols -1);
               
              // ... add central tet
              addTet(filtr, fcnvalues, curidx -1, curidx - ncols, curidx - nrows*ncols, curidx - nrows*ncols -ncols -1);
              // ... add remaining tets
              addOddTets(filtr, fcnvalues, ncols, nrows, i, j, k);
           }
        } // end for through j k positive
     } // end else through odd case.

    return;
} // end function addTriangles


int simplicesFromGrid(Fltr& filtr, const std::string& infile)
{
  std::ifstream in(infile.c_str());
  std::string	line;
  int nrows, ncols, ndimz;
  if (std::getline(in,line))
  {
    std::stringstream linestream(line);
    if (linestream >> nrows)
    {  if (linestream >> ncols)
       { if (linestream >> ndimz)
         {; //std::cout << nrows << " rows and " << ncols << " columns." << std::endl;
         }
          else
            ndimz = 1;  // to make backwards compatible with 2d grid
       }
    }
    else
      return 1;
  }
  
  int i = 0; // indexing the columns
  int j = 0; // indexing the rows
  int k = 0; // indexing the z dimension
  std::vector<double> fcnvalues;
  //double fcnvalues [i*j]; 

  while(std::getline(in, line))
  {
    //std::vector<Vertex> currow[ncols];
    if(line[0] == '#') continue;	// comment line
    std::stringstream	linestream(line);
    double x;
    while (linestream >> x) // each line corresponds to changing i
    {
      int curidx = i + ncols*j + ncols*nrows*k;
      fcnvalues.push_back(x); // at index i + ncols*j    
      assert(fcnvalues.at(curidx) == x);
		
      // .. add the vertex 
      std::vector<Vertex> vcont;
      vcont.push_back((Vertex)(curidx));
      filtr.push_back(Smplx(vcont, fcnvalues.at(curidx))); 

      // .. NEXT, Add the edges:
      addAllEdges(filtr, fcnvalues, ncols, nrows, i, j, k);
      // ... now add the triangles:
      addTriNTet(filtr, fcnvalues, ncols, nrows, i, j, k);

      ++i; // advance column
    } // end inner while loop, which iterates through i (a row / line)

    // ... advance row / z value
    i = 0;
    ++j;
    if (j > nrows -1)
    {
	j = 0;
        ++k;
    }
  } // end while 
  in.close();
  
  return 0;
} // end simplicesFromGrid function

#include <numeric>


int simplicesFromGridMem(Fltr & filtr, const double * const extFcnVal, const std::vector< unsigned int > & argGridNum)
{

    unsigned int idxCur, idxDim;
    const unsigned int gridNumProd = std::accumulate( argGridNum.begin(), argGridNum.end(), 1, std::multiplies< unsigned int >() );

	int ncols, nrows;
	ncols = nrows = 1;
  int i = 0; // indexing the columns
  int j = 0; // indexing the rows
  int k = 0; // indexing the z dimension
  std::vector<double> fcnvalues;
  int curidx = 0;

	if (argGridNum.size() > 0)
		ncols = argGridNum[0];
	if (argGridNum.size() > 1)
		nrows = argGridNum[1];

  while(curidx+1< gridNumProd )
  {
      curidx = i + argGridNum[0]*j + argGridNum[0]*argGridNum[1]*k;
      fcnvalues.push_back(extFcnVal[curidx]); // at index i + ncols*j    

      // .. add the vertex 
      std::vector<Vertex> vcont;
      vcont.push_back((Vertex)(curidx));
      filtr.push_back(Smplx(vcont, fcnvalues.at(curidx))); 

      // .. NEXT, Add the edges:
      addAllEdges(filtr, fcnvalues, ncols, nrows, i, j, k);
      // ... now add the triangles:
      addTriNTet(filtr, fcnvalues, ncols, nrows, i, j, k);

      ++i; // advance column

     // ... advance row / z value
     if (i >= argGridNum[0])
     {
    	i = 0;
    	++j;
     }
     if (j >= argGridNum[1])
     {
	j = 0;
        ++k;
     }

  }

  
  return 0;
} // end simplicesFromGrid function

inline std::vector< unsigned int > gridNumber(const int * const extDim, const int * const extGridNum)
{
    return std::vector< unsigned int >( extGridNum, extGridNum+extDim[0] );
}
