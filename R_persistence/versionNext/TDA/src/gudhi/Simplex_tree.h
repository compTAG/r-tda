/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Clément Maria
 *
 *    Copyright (C) 2014  INRIA Sophia Antipolis-Méditerranée (France)
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef SRC_SIMPLEX_TREE_INCLUDE_GUDHI_SIMPLEX_TREE_H_
#define SRC_SIMPLEX_TREE_INCLUDE_GUDHI_SIMPLEX_TREE_H_

#include <gudhi/Simplex_tree/Simplex_tree_node_explicit_storage.h>
#include <gudhi/Simplex_tree/Simplex_tree_siblings.h>
#include <gudhi/Simplex_tree/Simplex_tree_iterators.h>
#include <gudhi/Simplex_tree/indexing_tag.h>

#include <boost/container/flat_map.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <boost/graph/adjacency_list.hpp>

#include <algorithm>
#include <utility>
#include <vector>
#include <limits>


namespace Gudhi {
  
/** \defgroup simplex_tree Filtered Complexes
/** \defgroup simplex_tree Filtered Complexes
 *
 * A simplicial complex \f$\mathbf{K}\f$
 * on a set of vertices \f$V = \{1, \cdots ,|V|\}\f$ is a collection of simplices
 * \f$\{\sigma\}\f$,
 * \f$\sigma \subseteq V\f$ such that \f$\tau \subseteq \sigma \in \mathbf{K} \rightarrow \tau \in
 * \mathbf{K}\f$. The
 * dimension \f$n=|\sigma|-1\f$ of \f$\sigma\f$ is its number of elements minus \f$1\f$.
 *
 * A filtration of a simplicial complex is
 * a function \f$f:\mathbf{K} \rightarrow \mathbb{R}\f$ satisfying \f$f(\tau)\leq f(\sigma)\f$ whenever
 * \f$\tau \subseteq \sigma\f$. Ordering the simplices by increasing filtration values
 * (breaking ties so as a simplex appears after its subsimplices of same filtration value)
 * provides an indexing scheme.
 *

 <DT>Implementations:</DT>
 There are two implementation of complexes. The first on is the Simplex_tree data structure.
 The simplex tree is an efficient and flexible
 data structure for representing general (filtered) simplicial complexes. The data structure
 is described in \cite boissonnatmariasimplextreealgorithmica

 The second one is the Hasse_complex. The Hasse complex is a data structure representing
 explicitly all co-dimension 1 incidence relations in a complex. It is consequently faster
 when accessing the boundary of a simplex, but is less compact and harder to construct from
 scratch.


 * \author    Clément Maria
 * \version   1.0
 * \date      2014
 * \copyright GNU General Public License v3.
 * @{
 */
/**
 * \brief Simplex Tree data structure for representing simplicial complexes.
 *
 * \details Every simplex \f$[v_0, \cdots ,v_d]\f$ admits a canonical orientation
 * induced by the order relation on vertices \f$ v_0 < \cdots < v_d \f$.
 *
 * Details may be found in \cite boissonnatmariasimplextreealgorithmica.
 *
 * \implements FilteredComplex
 *
 */
template<typename IndexingTag = linear_indexing_tag,
    typename FiltrationValue = double, typename SimplexKey = int  // must be a signed integer type
    , typename VertexHandle = int  // must be a signed integer type, int convertible to it
//         , bool ContiguousVertexHandles = true   //true is Vertex_handles are exactly the set [0;n)
>
class Simplex_tree {
 public:
  typedef IndexingTag Indexing_tag;
  /** \brief Type for the value of the filtration function.
   *
   * Must be comparable with <. */
  typedef FiltrationValue Filtration_value;
  /** \brief Key associated to each simplex.
   *
   * Must be a signed integer type. */
  typedef SimplexKey Simplex_key;
  /** \brief Type for the vertex handle.
   *
   * Must be a signed integer type. It admits a total order <. */
  typedef VertexHandle Vertex_handle;

  /* Type of node in the simplex tree. */
  typedef Simplex_tree_node_explicit_storage<Simplex_tree> Node;
  /* Type of dictionary Vertex_handle -> Node for traversing the simplex tree. */
  typedef typename boost::container::flat_map<Vertex_handle, Node> Dictionary;

  friend class Simplex_tree_node_explicit_storage< Simplex_tree<FiltrationValue, SimplexKey, VertexHandle> >;
  friend class Simplex_tree_siblings< Simplex_tree<FiltrationValue, SimplexKey, VertexHandle>, Dictionary>;
  friend class Simplex_tree_simplex_vertex_iterator< Simplex_tree<FiltrationValue, SimplexKey, VertexHandle> >;
  friend class Simplex_tree_boundary_simplex_iterator< Simplex_tree<FiltrationValue, SimplexKey, VertexHandle> >;
  friend class Simplex_tree_complex_simplex_iterator< Simplex_tree<FiltrationValue, SimplexKey, VertexHandle> >;
  friend class Simplex_tree_skeleton_simplex_iterator< Simplex_tree<FiltrationValue, SimplexKey, VertexHandle> >;

  /* \brief Set of nodes sharing a same parent in the simplex tree. */
  /* \brief Set of nodes sharing a same parent in the simplex tree. */
  typedef Simplex_tree_siblings<Simplex_tree, Dictionary> Siblings;

 public:
  /** \brief Handle type to a simplex contained in the simplicial complex represented
   * byt he simplex tree. */
  typedef typename Dictionary::iterator Simplex_handle;

 private:
  typedef typename Dictionary::iterator Dictionary_it;
  typedef typename Dictionary_it::value_type Dit_value_t;

  struct return_first {
    Vertex_handle operator()(const Dit_value_t& p_sh) const {
      return p_sh.first;
    }
  };

 public:
  /** \name Range and iterator types
   *
   * The naming convention is Container_content_(iterator/range). A Container_content_range is
   * essentially an object on which the methods begin() and end() can be called. They both return
   * an object of type Container_content_iterator, and allow the traversal of the range
   * [ begin();end() ).
   * @{ */

  /** \brief Iterator over the vertices of the simplicial complex.
   *
   * 'value_type' is Vertex_handle. */
  typedef boost::transform_iterator<return_first, Dictionary_it> Complex_vertex_iterator;
  /** \brief Range over the vertices of the simplicial complex. */
  typedef boost::iterator_range<Complex_vertex_iterator> Complex_vertex_range;
  /** \brief Iterator over the vertices of a simplex.
   *
   * 'value_type' is Vertex_handle. */
  typedef Simplex_tree_simplex_vertex_iterator<Simplex_tree> Simplex_vertex_iterator;
  /** \brief Range over the vertices of a simplex. */
  typedef boost::iterator_range<Simplex_vertex_iterator> Simplex_vertex_range;
  /** \brief Iterator over the simplices of the boundary of a simplex.
   *
   * 'value_type' is Simplex_handle. */
  typedef Simplex_tree_boundary_simplex_iterator<Simplex_tree> Boundary_simplex_iterator;
  /** \brief Range over the simplices of the boundary of a simplex. */
  typedef boost::iterator_range<Boundary_simplex_iterator> Boundary_simplex_range;
  /** \brief Iterator over the simplices of the simplicial complex.
   *
   * 'value_type' is Simplex_handle. */
  typedef Simplex_tree_complex_simplex_iterator<Simplex_tree> Complex_simplex_iterator;
  /** \brief Range over the simplices of the simplicial complex. */
  typedef boost::iterator_range<Complex_simplex_iterator> Complex_simplex_range;
  /** \brief Iterator over the simplices of the skeleton of the simplicial complex, for a given 
   * dimension.
   *
   * 'value_type' is Simplex_handle. */
  typedef Simplex_tree_skeleton_simplex_iterator<Simplex_tree> Skeleton_simplex_iterator;
  /** \brief Range over the simplices of the skeleton of the simplicial complex, for a given 
   * dimension. */
  typedef boost::iterator_range<Skeleton_simplex_iterator> Skeleton_simplex_range;
  /** \brief Iterator over the simplices of the simplicial complex, ordered by the filtration.
   *
   * 'value_type' is Simplex_handle. */
  typedef typename std::vector<Simplex_handle>::iterator Filtration_simplex_iterator;
  /** \brief Range over the simplices of the simplicial complex, ordered by the filtration. */
  typedef boost::iterator_range<Filtration_simplex_iterator> Filtration_simplex_range;

  /* @} */  // end name range and iterator types
  /** \name Range and iterator methods
   * @{ */

  /** \brief Returns a range over the vertices of the simplicial complex.
   *
   * The order is increasing according to < on Vertex_handles.*/
  Complex_vertex_range complex_vertex_range() {
    return Complex_vertex_range(
        boost::make_transform_iterator(root_.members_.begin(), return_first()),
        boost::make_transform_iterator(root_.members_.end(), return_first()));
  }

  /** \brief Returns a range over the simplices of the simplicial complex.
   *
   * In the Simplex_tree, the tree is traverse in a depth-first fashion.
   * Consequently, simplices are ordered according to lexicographic order on the list of
   * Vertex_handles of a simplex, read in increasing < order for Vertex_handles. */
  Complex_simplex_range complex_simplex_range() {
    return Complex_simplex_range(Complex_simplex_iterator(this),
                                 Complex_simplex_iterator());
  }

  /** \brief Returns a range over the simplices of the dim-skeleton of the simplicial complex.
   *
   * The \f$d\f$-skeleton of a simplicial complex \f$\mathbf{K}\f$ is the simplicial complex containing the
   * simplices of \f$\mathbf{K}\f$ of dimension at most \f$d\f$.
   *
   * @param[in] dim The maximal dimension of the simplices in the skeleton.
   *
   * The simplices are ordered according to lexicographic order on the list of
   * Vertex_handles of a simplex, read in increasing < order for Vertex_handles. */
  Skeleton_simplex_range skeleton_simplex_range(int dim) {
    return Skeleton_simplex_range(Skeleton_simplex_iterator(this, dim),
                                  Skeleton_simplex_iterator());
  }

  /** \brief Returns a range over the simplices of the simplicial complex,
   * in the order of the filtration.
   *
   * The filtration is a monotonic function \f$ f: \mathbf{K} \rightarrow \mathbb{R} \f$, i.e. if two simplices
   * \f$\tau\f$ and \f$\sigma\f$ satisfy \f$\tau \subseteq \sigma\f$ then
   * \f$f(\tau) \leq f(\sigma)\f$.
   *
   * The method returns simplices ordered according to increasing filtration values. Ties are
   * resolved by considering inclusion relation (subsimplices appear before their cofaces). If two
   * simplices have same filtration value but are not comparable w.r.t. inclusion, lexicographic
   * order is used.
   *
   * The filtration must be valid. If the filtration has not been initialized yet, the
   * method initializes it (i.e. order the simplices). */
  Filtration_simplex_range filtration_simplex_range(linear_indexing_tag) {
    if (filtration_vect_.empty()) {
      initialize_filtration();
    }
    return Filtration_simplex_range(filtration_vect_.begin(),
                                    filtration_vect_.end());
  }

  Filtration_simplex_range filtration_simplex_range() {
    return filtration_simplex_range(Indexing_tag());
  }
  /** \brief Returns a range over the vertices of a simplex.
   *
   * The order in which the vertices are visited is the decreasing order for < on Vertex_handles,
   * which is consequenlty
   * equal to \f$(-1)^{\text{dim} \sigma}\f$ the canonical orientation on the simplex.
   */
  Simplex_vertex_range simplex_vertex_range(Simplex_handle sh) {
    return Simplex_vertex_range(Simplex_vertex_iterator(this, sh),
                                Simplex_vertex_iterator(this));
  }

  /** \brief Returns a range over the simplices of the boundary of a simplex.
   *
   * The boundary of a simplex is the set of codimension \f$1\f$ subsimplices of the simplex.
   * If the simplex is \f$[v_0, \cdots ,v_d]\f$, with canonical orientation
   * induced by \f$ v_0 < \cdots < v_d \f$, the iterator enumerates the
   * simplices of the boundary in the order:
   * \f$[v_0,\cdots,\widehat{v_i},\cdots,v_d]\f$ for \f$i\f$ from \f$0\f$ to \f$d\f$,
   * where \f$\widehat{v_i}\f$ means that the vertex \f$v_i\f$ is omitted.
   *
   * We note that the alternate sum of the simplices given by the iterator
   * gives \f$(-1)^{\text{dim} \sigma}\f$ the chains corresponding to the boundary
   * of the simplex.
   *
   * @param[in] sh Simplex for which the boundary is computed. */
  Boundary_simplex_range boundary_simplex_range(Simplex_handle sh) {
    return Boundary_simplex_range(Boundary_simplex_iterator(this, sh),
                                  Boundary_simplex_iterator(this));
  }

  /** @} */  // end range and iterator methods
  /** \name Constructor/Destructor
   * @{ */

  /** \brief Constructs an empty simplex tree. */
  Simplex_tree()
      : null_vertex_(-1),
        threshold_(0),
        num_simplices_(0),
        root_(NULL, null_vertex_),
        filtration_vect_(),
        dimension_(-1) {
  }

  /** \brief Destructor; deallocates the whole tree structure. */
  ~Simplex_tree() {
    for (auto sh = root_.members().begin(); sh != root_.members().end(); ++sh) {
      if (has_children(sh)) {
        rec_delete(sh->second.children());
      }
    }
  }
  /** @} */  // end constructor/destructor
 private:
  /** Recursive deletion. */
  void rec_delete(Siblings * sib) {
    for (auto sh = sib->members().begin(); sh != sib->members().end(); ++sh) {
      if (has_children(sh)) {
        rec_delete(sh->second.children());
      }
    }
    delete sib;
  }

 public:
  /** \brief Returns the key associated to a simplex.
   *
   * The filtration must be initialized. */
  Simplex_key key(Simplex_handle sh) {
    return sh->second.key();
  }
  /** \brief Returns the simplex associated to a key.
   *
   * The filtration must be initialized. */
  Simplex_handle simplex(Simplex_key key) {
    return filtration_vect_[key];
  }
  /** \brief Returns the filtration value of a simplex.
   *
   * Called on the null_simplex, returns std::numeric_limits<double>::infinity(). 
   * std::numeric_limits<Filtration_value>::has_infinity must be 1.*/
  Filtration_value filtration(Simplex_handle sh) { 
    if(sh != null_simplex()) { return sh->second.filtration(); }
    else                     { return std::numeric_limits< Filtration_value >::infinity(); }//filtration(); }
  }
  /** \brief Returns an upper bound of the filtration values of the simplices. */
  Filtration_value filtration() {
    return threshold_;
  }
  /** \brief Returns a Simplex_handle different from all Simplex_handles
   * associated to the simplices in the simplicial complex.
   *
   * One can call filtration(null_simplex()). */
  Simplex_handle null_simplex() {
    return Dictionary_it(NULL);
  }
  /** \brief Returns a key different for all keys associated to the
   * simplices of the simplicial complex. */
  Simplex_key null_key() {
    return -1;
  }
  /** \brief Returns a Vertex_handle different from all Vertex_handles associated
   * to the vertices of the simplicial complex. */
  Vertex_handle null_vertex() {
    return null_vertex_;
  }
  /** \brief Returns the number of vertices in the complex. */
  size_t num_vertices() {
    return root_.members_.size();
  }
  /** \brief Returns the number of simplices in the complex.
   *
   * Does not count the empty simplex. */
  const unsigned int& num_simplices() const {
    return num_simplices_;
  }

  /** \brief Returns the dimension of a simplex.
   *
   * Must be different from null_simplex().*/
  int dimension(Simplex_handle sh) {
    Siblings * curr_sib = self_siblings(sh);
    int dim = 0;
    while (curr_sib != NULL) {
      ++dim;
      curr_sib = curr_sib->oncles();
    }
    return dim - 1;
  }
  /** \brief Returns an upper bound on the dimension of the simplicial complex. */
  int dimension() {
    return dimension_;
  }

  /** \brief Returns true iff the node in the simplex tree pointed by
   * sh has children.*/
  bool has_children(Simplex_handle sh) {
    return (sh->second.children()->parent() == sh->first);
  }

  /** \brief Given a range of Vertex_handles, returns the Simplex_handle
   * of the simplex in the simplicial complex containing the corresponding
   * vertices. Return null_simplex() if the simplex is not in the complex.
   *
   * The type RandomAccessVertexRange must be a range for which .begin() and
   * .end() return random access iterators, with <CODE>value_type</CODE>
   * <CODE>Vertex_handle</CODE>.
   */
  template<class RandomAccessVertexRange>
  Simplex_handle find(const RandomAccessVertexRange & s) {
#ifdef DEBUG_TRACES
    if (s.begin() == s.end())
      std::cerr << "Empty simplex \n";
#endif

    sort(s.begin(), s.end());

    Siblings * tmp_sib = &root_;
    Dictionary_it tmp_dit;
    Vertex_handle last = s[s.size() - 1];
    for (auto v : s) {
      tmp_dit = tmp_sib->members_.find(v);
      if (tmp_dit == tmp_sib->members_.end()) {
        return null_simplex();
      }
      if (!has_children(tmp_dit) && v != last) {
        return null_simplex();
      }
      tmp_sib = tmp_dit->second.children();
    }
    return tmp_dit;
  }

  /** \brief Returns the Simplex_handle corresponding to the 0-simplex
   * representing the vertex with Vertex_handle v. */
  Simplex_handle find_vertex(Vertex_handle v) {
    return root_.members_.begin() + v;
  }
//{ return root_.members_.find(v); }

  /** \brief Insert a simplex, represented by a range of Vertex_handles, in the simplicial complex.
   *
   * @param[in]  simplex    range of Vertex_handles, representing the vertices of the new simplex
   * @param[in]  filtration the filtration value assigned to the new simplex.
   * The return type is a pair. If the new simplex is inserted successfully (i.e. it was not in the
   * simplicial complex yet) the bool is set to true and the Simplex_handle is the handle assigned
   * to the new simplex.
   * If the insertion fails (the simplex is already there), the bool is set to false. If the insertion
   * fails and the simplex already in the complex has a filtration value strictly bigger than 'filtration',
   * we assign this simplex with the new value 'filtration', and set the Simplex_handle filed of the
   * output pair to the Simplex_handle of the simplex. Otherwise, we set the Simplex_handle part to
   * null_simplex.
   *
   * All subsimplices do not necessary need to be already in the simplex tree to proceed to an
   * insertion. However, the property of being a simplicial complex will be violated. This allows
   * us to insert a stream of simplices contained in a simplicial complex without considering any
   * order on them.
   *
   * The filtration value
   * assigned to the new simplex must preserve the monotonicity of the filtration.
   *
   * The type RandomAccessVertexRange must be a range for which .begin() and
   * .end() return random access iterators, with 'value_type' Vertex_handle. */
  template<class RandomAccessVertexRange>
  std::pair<Simplex_handle, bool> insert(RandomAccessVertexRange & simplex,
                                         Filtration_value filtration) {
    if (simplex.empty()) {
      return std::pair<Simplex_handle, bool>(null_simplex(), true);
    }

    sort(simplex.begin(), simplex.end());  // must be sorted in increasing order

    Siblings * curr_sib = &root_;
    std::pair<Simplex_handle, bool> res_insert;
    typename RandomAccessVertexRange::iterator vi;
    for (vi = simplex.begin(); vi != simplex.end() - 1; ++vi) {
      res_insert = curr_sib->members_.emplace(*vi, Node(curr_sib, filtration));
      if (!(has_children(res_insert.first))) {
        res_insert.first->second.assign_children(new Siblings(curr_sib, *vi));
      }
      curr_sib = res_insert.first->second.children();
    }
    res_insert = curr_sib->members_.emplace(*vi, Node(curr_sib, filtration));
    if (!res_insert.second) {  // if already in the complex
      if (res_insert.first->second.filtration() > filtration) {  // if filtration value modified
        res_insert.first->second.assign_filtration(filtration);
        return res_insert;
      }
      return std::pair<Simplex_handle, bool>(null_simplex(), false);  // if filtration value unchanged
    }
    // otherwise the insertion has succeeded
    return res_insert;
  }

  /** \brief Assign a value 'key' to the key of the simplex
   * represented by the Simplex_handle 'sh'. */
  void assign_key(Simplex_handle sh, Simplex_key key) {
    sh->second.assign_key(key);
  }

 public:
  /** Returns the two Simplex_handle corresponding to the endpoints of
   * and edge. sh must point to a 1-dimensional simplex. This is an
   * optimized version of the boundary computation. */
  std::pair<Simplex_handle, Simplex_handle> endpoints(Simplex_handle sh) {
    return std::pair<Simplex_handle, Simplex_handle>(
        root_.members_.find(sh->first),
        root_.members_.find(self_siblings(sh)->parent()));
  }

  /** Returns the Siblings containing a simplex.*/
  Siblings * self_siblings(Simplex_handle sh) {
    if (sh->second.children()->parent() == sh->first)
      return sh->second.children()->oncles();
    else
      return sh->second.children();
  }

// void display_simplex(Simplex_handle sh)
// {
//   std::cout << "   " << "[" << filtration(sh) << "] ";
//   for( auto vertex : simplex_vertex_range(sh) )
//     { std::cout << vertex << " "; }
// }

  // void print(Simplex_handle sh, std::ostream& os = std::cout)
  // {  for(auto v : simplex_vertex_range(sh)) {os << v << " ";}
  //    os << std::endl; }

 public:
  /** Returns a pointer to the root nodes of the simplex tree. */
  Siblings * root() {
    return &root_;
  }

 public:
  /** Set an upper bound for the filtration values. */
  void set_filtration(Filtration_value fil) {
    threshold_ = fil;
  }
  /** Set a number of simplices for the simplicial complex. */
  void set_num_simplices(const unsigned int& num_simplices) {
    num_simplices_ = num_simplices;
  }
  /** Set a dimension for the simplicial complex. */
  void set_dimension(int dimension) {
    dimension_ = dimension;
  }

 public:
  /** \brief Initializes the filtrations, i.e. sort the
   * simplices according to their order in the filtration and initializes all Simplex_keys.
   *
   * After calling this method, filtration_simplex_range() becomes valid, and each simplex is
   * assigned a Simplex_key corresponding to its order in the filtration (from 0 to m-1 for a
   * simplicial complex with m simplices).
   *
   * The use of a depth-first traversal of the simplex tree, provided by
   * complex_simplex_range(), combined with
   * a stable sort is meant to optimize the order of simplices with same
   * filtration value. The heuristic consists in inserting the cofaces of a
   * simplex as soon as possible.
   *
   * Will be automatically called when calling filtration_simplex_range()
   * if the filtration has never been initialized yet. */
  void initialize_filtration() {
    filtration_vect_.clear();
    filtration_vect_.reserve(num_simplices());
    for (auto cpx_it = complex_simplex_range().begin();
        cpx_it != complex_simplex_range().end(); ++cpx_it) {
      filtration_vect_.push_back(*cpx_it);
    }

    // the stable sort ensures the ordering heuristic
    std::stable_sort(filtration_vect_.begin(), filtration_vect_.end(),
                     is_before_in_filtration(this));
  }

 private:
  /** \brief Returns true iff the list of vertices of sh1
   * is smaller than the list of vertices of sh2 w.r.t.
   * lexicographic order on the lists read in reverse.
   *
   * It defines a StrictWeakOrdering on simplices. The Simplex_vertex_iterators
   * must traverse the Vertex_handle in decreasing order. Reverse lexicographic order satisfy
   * the property that a subsimplex of a simplex is always strictly smaller with this order. */
  bool reverse_lexicographic_order(Simplex_handle sh1, Simplex_handle sh2) {
    Simplex_vertex_range rg1 = simplex_vertex_range(sh1);
    Simplex_vertex_range rg2 = simplex_vertex_range(sh2);
    Simplex_vertex_iterator it1 = rg1.begin();
    Simplex_vertex_iterator it2 = rg2.begin();
    while (it1 != rg1.end() && it2 != rg2.end()) {
      if (*it1 == *it2) {
        ++it1;
        ++it2;
      } else {
        return *it1 < *it2;
      }
    }
    return ((it1 == rg1.end()) && (it2 != rg2.end()));
  }
  /** \brief StrictWeakOrdering, for the simplices, defined by the filtration.
   *
   * It corresponds to the partial order
   * induced by the filtration values, with ties resolved using reverse lexicographic order.
   * Reverse lexicographic order has the property to always consider the subsimplex of a simplex
   * to be smaller. The filtration function must be monotonic. */
  struct is_before_in_filtration {
    explicit is_before_in_filtration(Simplex_tree * st)
        : st_(st) {
    }

    bool operator()(const Simplex_handle sh1, const Simplex_handle sh2) const {
      if (st_->filtration(sh1) != st_->filtration(sh2)) {
        return st_->filtration(sh1) < st_->filtration(sh2);
      }

      return st_->reverse_lexicographic_order(sh1, sh2);  // is sh1 a proper subface of sh2
    }

    Simplex_tree * st_;
  };

 public:
  /** \brief Inserts a 1-skeleton in an empty Simplex_tree.
   *
   * The Simplex_tree must contain no simplex when the method is
   * called.
   *
   * Inserts all vertices and edges given by a OneSkeletonGraph.
   * OneSkeletonGraph must be a model of boost::AdjacencyGraph,
   * boost::EdgeListGraph and boost::PropertyGraph.
   *
   * The vertex filtration value is accessible through the property tag
   * vertex_filtration_t.
   * The edge filtration value is accessible through the property tag
   * edge_filtration_t.
   *
   * boost::graph_traits<OneSkeletonGraph>::vertex_descriptor
   *                                    must be Vertex_handle.
   * boost::graph_traits<OneSkeletonGraph>::directed_category
   *                                    must be undirected_tag. */
  template<class OneSkeletonGraph>
  void insert_graph(const OneSkeletonGraph& skel_graph) {
    assert(num_simplices() == 0);  // the simplex tree must be empty

    if (boost::num_vertices(skel_graph) == 0) {
      return;
    }
    if (num_edges(skel_graph) == 0) {
      dimension_ = 0;
    } else {
      dimension_ = 1;
    }

    num_simplices_ = boost::num_vertices(skel_graph)
        + boost::num_edges(skel_graph);
    root_.members_.reserve(boost::num_vertices(skel_graph));

    typename boost::graph_traits<OneSkeletonGraph>::vertex_iterator v_it,
        v_it_end;
    for (std::tie(v_it, v_it_end) = boost::vertices(skel_graph); v_it != v_it_end;
        ++v_it) {
      root_.members_.emplace_hint(
          root_.members_.end(), *v_it,
          Node(&root_, boost::get(vertex_filtration_t(), skel_graph, *v_it)));
    }
    typename boost::graph_traits<OneSkeletonGraph>::edge_iterator e_it,
        e_it_end;
    for (std::tie(e_it, e_it_end) = boost::edges(skel_graph); e_it != e_it_end;
        ++e_it) {
      auto u = source(*e_it, skel_graph);
      auto v = target(*e_it, skel_graph);
      if (u < v) {  // count edges only once { std::swap(u,v); } // u < v
        auto sh = find_vertex(u);
        if (!has_children(sh)) {
          sh->second.assign_children(new Siblings(&root_, sh->first));
        }

        sh->second.children()->members().emplace(
            v,
            Node(sh->second.children(),
                 boost::get(edge_filtration_t(), skel_graph, *e_it)));
      }
    }
  }
  /** \brief Expands the Simplex_tree containing only its one skeleton
   * until dimension max_dim.
   *
   * The expanded simplicial complex until dimension \f$d\f$
   * attached to a graph \f$G\f$ is the maximal simplicial complex of
   * dimension at most \f$d\f$ admitting the graph \f$G\f$ as \f$1\f$-skeleton.
   * The filtration value assigned to a simplex is the maximal filtration
   * value of one of its edges.
   *
   * The Simplex_tree must contain no simplex of dimension bigger than
   * 1 when calling the method. */
  void expansion(int max_dim) {
    dimension_ = max_dim;
    for (Dictionary_it root_it = root_.members_.begin();
        root_it != root_.members_.end(); ++root_it) {
      if (has_children(root_it)) {
        siblings_expansion(root_it->second.children(), max_dim - 1);
      }
    }
    dimension_ = max_dim - dimension_;
  }

 private:
  /** \brief Recursive expansion of the simplex tree.*/
  void siblings_expansion(Siblings * siblings,  // must contain elements
      int k) {
    if (dimension_ > k) {
      dimension_ = k;
    }
    if (k == 0)
      return;
    Dictionary_it next = siblings->members().begin();
    ++next;

    static std::vector<std::pair<Vertex_handle, Node> > inter;  // static, not thread-safe.
    for (Dictionary_it s_h = siblings->members().begin();
        s_h != siblings->members().end(); ++s_h, ++next) {
      Simplex_handle root_sh = find_vertex(s_h->first);
      if (has_children(root_sh)) {
        intersection(
            inter,                      // output intersection
            next,                       // begin
            siblings->members().end(),  // end
            root_sh->second.children()->members().begin(),
            root_sh->second.children()->members().end(),
            s_h->second.filtration());
        if (inter.size() != 0) {
          this->num_simplices_ += inter.size();
          Siblings * new_sib = new Siblings(siblings,  // oncles
              s_h->first,  // parent
              inter);      // boost::container::ordered_unique_range_t
          inter.clear();
          s_h->second.assign_children(new_sib);
          siblings_expansion(new_sib, k - 1);
        } else {
          s_h->second.assign_children(siblings);  // ensure the children property
          inter.clear();
        }
      }
    }
  }
  /** \brief Intersects Dictionary 1 [begin1;end1) with Dictionary 2 [begin2,end2)
   * and assigns the maximal possible Filtration_value to the Nodes. */
  void intersection(std::vector<std::pair<Vertex_handle, Node> >& intersection,
                    Dictionary_it begin1, Dictionary_it end1,
                    Dictionary_it begin2, Dictionary_it end2,
                    Filtration_value filtration) {
    if (begin1 == end1 || begin2 == end2)
      return;  // ----->>
    while (true) {
      if (begin1->first == begin2->first) {
        intersection.push_back(
            std::pair<Vertex_handle, Node>(
                begin1->first,
                Node(NULL, maximum(begin1->second.filtration(), begin2->second.filtration(), filtration))));
        ++begin1;
        ++begin2;
        if (begin1 == end1 || begin2 == end2)
          return;  // ----->>
      } else {
        if (begin1->first < begin2->first) {
          ++begin1;
          if (begin1 == end1)
            return;
        } else {
          ++begin2;
          if (begin2 == end2)
            return;  // ----->>
        }
      }
    }
  }
  /** Maximum over 3 values.*/
  Filtration_value maximum(Filtration_value a, Filtration_value b,
                           Filtration_value c) {
    Filtration_value max = (a < b) ? b : a;
    return ((max < c) ? c : max);
  }

 public:
  /** \brief Write the hasse diagram of the simplicial complex in os.
   *
   * Each row in the file correspond to a simplex. A line is written:
   * dim idx_1 ... idx_k fil   where dim is the dimension of the simplex,
   * idx_1 ... idx_k are the row index (starting from 0) of the simplices of the boundary
   * of the simplex, and fil is its filtration value. */
  void print_hasse(std::ostream& os) {
    os << num_simplices() << " " << std::endl;
    for (auto sh : filtration_simplex_range(Indexing_tag())) {
      os << dimension(sh) << " ";
      for (auto b_sh : boundary_simplex_range(sh)) {
        os << key(b_sh) << " ";
      }
      os << filtration(sh) << " \n";
    }
  }

 private:
  Vertex_handle null_vertex_;
  /** \brief Upper bound on the filtration values of the simplices.*/
  Filtration_value threshold_;
  /** \brief Total number of simplices in the complex, without the empty simplex.*/
  unsigned int num_simplices_;
  /** \brief Set of simplex tree Nodes representing the vertices.*/
  Siblings root_;
  /** \brief Simplices ordered according to a filtration.*/
  std::vector<Simplex_handle> filtration_vect_;
  /** \brief Upper bound on the dimension of the simplicial complex.*/
  int dimension_;
};

// Print a Simplex_tree in os.
template<typename T1, typename T2, typename T3>
std::ostream& operator<<(std::ostream & os, Simplex_tree<T1, T2, T3> & st) {
  for (auto sh : st.filtration_simplex_range()) {
    os << st.dimension(sh) << " ";
    for (auto v : st.simplex_vertex_range(sh)) {
      os << v << " ";
    }
    os << st.filtration(sh) << "\n";  // TODO(VR): why adding the key ?? not read ?? << "     " << st.key(sh) << " \n";
  }
  return os;
}
template<typename T1, typename T2, typename T3>
std::istream& operator>>(std::istream & is, Simplex_tree<T1, T2, T3> & st) {
  // assert(st.num_simplices() == 0);

  std::vector<typename Simplex_tree<T1, T2, T3>::Vertex_handle> simplex;
  typename Simplex_tree<T1, T2, T3>::Filtration_value fil;
  typename Simplex_tree<T1, T2, T3>::Filtration_value max_fil = 0;
  int max_dim = -1;
  size_t num_simplices = 0;
  while (read_simplex(is, simplex, fil)) {  // read all simplices in the file as a list of vertices
    ++num_simplices;
    int dim = static_cast<int>(simplex.size() - 1);  // Warning : simplex_size needs to be casted in int - Can be 0
    if (max_dim < dim) {
      max_dim = dim;
    }
    if (max_fil < fil) {
      max_fil = fil;
    }
    st.insert(simplex, fil);  // insert every simplex in the simplex tree
    simplex.clear();
  }
  st.set_num_simplices(num_simplices);
  st.set_dimension(max_dim);
  st.set_filtration(max_fil);

  return is;
}

/** @} */  // end defgroup simplex_tree
}  // namespace Gudhi

#endif  // SRC_SIMPLEX_TREE_INCLUDE_GUDHI_SIMPLEX_TREE_H_
