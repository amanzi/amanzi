#ifndef _COORDINATES_HH_
#define _COORDINATES_HH_

#include "dbc.hh"
#include <vector>
#include <memory>

namespace Mesh_data
{

template <typename S>
class Coordinates
{

    int num_nodes_;
    int dimension_;

    std::vector< std::vector<S> > coordinates_;

    void init_storage_ ();

    template <typename IT>
    void assign_coord_ (IT it, int dimension);

    template <typename IT>
    void assign_data_ (IT it_x, IT it_y, IT it_z);

    bool ok_dim_  (int dim)  const { return dim >= 0  && dim < dimension_; }
    bool ok_node_ (int node) const { return node >= 0 && node < num_nodes_; }

    Coordinates<S>(const Coordinates<S>& rhs);

    bool valid () const;

public:

    // Information
    // -----------

    typedef S value_type;
    typedef S* pointer_type;
    typedef S& reference_type;

    int nodes ()     const { return num_nodes_; }
    int dimension () const { return dimension_; }


    // Constructors & Operators
    // ------------------------

    Coordinates () : num_nodes_ (0), dimension_ (0), coordinates_ () { }
    Coordinates (int num_nodes, int dimension) : num_nodes_ (num_nodes), dimension_ (dimension), 
                                                 coordinates_ (dimension) 
    { init_storage_ ();  }

    template <typename IT>
    Coordinates (int num_nodes, IT iterator_x, IT iterator_y = IT (0), IT iterator_z = IT (0));

    void swap (Coordinates<S>& rhs);

    template <typename T>
    void take_data_from (std::vector< std::vector<T> >& rhs);
    void take_data_from (std::vector< std::vector<S> >& rhs);

    template <typename T>
    static Coordinates<S>* build_from (std::vector< std::vector<T> > &rhs);

    void to_stream (std::ostream& stream, bool verbose = false) const;


    // Data accessors
    // --------------

    S& operator () (int node, int dimension);
    S  operator () (int node, int dimension) const;

    template <typename IT>
    void operator () (int node, IT store) const;

    S* get_coordinate_pointer (int dimension);
    S const * get_coordinate_pointer (int dimension) const;


};

// -----------------
// Private functions
// -----------------


template <typename S>
void Coordinates<S>::init_storage_ ()
{
    coordinates_.resize (dimension_);

    for (int d=0; d<dimension_; ++d)
    {
        coordinates_[d].resize(num_nodes_);
    }
}

template <typename S>
template <typename IT>
void Coordinates<S>::assign_coord_ (IT it, int dimension)
{
    if (it != IT (0))
    {
        ASSERT (ok_dim_ (dimension));
        coordinates_ [dimension].assign (it, it+num_nodes_);
    }
}

template <typename S>
template <typename IT>
void Coordinates<S>::assign_data_ (IT it_x, IT it_y, IT it_z)
{
    assign_coord_ (it_x, 0);
    assign_coord_ (it_y, 1);
    assign_coord_ (it_z, 2);
}


// -------------
// Data accessors
// -------------


template <typename S>
S& Coordinates<S>::operator () (int node, int dimension)
{
    ASSERT (ok_dim_ (dimension));
    ASSERT (ok_node_ (node));

    return coordinates_ [dimension] [node];
}

template <typename S>
S Coordinates<S>::operator () (int node, int dimension) const
{
    ASSERT (ok_dim_ (dimension));
    ASSERT (ok_node_ (node));

    return coordinates_ [dimension] [node];
}

template <typename S>
template <typename IT>
void Coordinates<S>::operator () (int node, IT coords) const
{
    ASSERT (ok_node_ (node));

    for (int d=0; d<dimension_; ++d)
    {
        *coords = (*this) (node, d);
        ++coords;
    }

}

template <typename S>
S* Coordinates<S>::get_coordinate_pointer (int dimension)
{
    if (ok_dim_ (dimension))
        return &( coordinates_ [dimension]) [0];

    return (S*)(0);
}

template <typename S>
S const * Coordinates<S>::get_coordinate_pointer (int dimension) const
{
    if (ok_dim_ (dimension))
        return &(coordinates_ [dimension]) [0];

    return (S const*) (0);
}



// ------------------------
// Constructors & Operators
// ------------------------

template <typename S>
template <typename IT>
Coordinates<S>::Coordinates (int num_nodes, IT iterator_x, IT iterator_y, IT iterator_z) : 
    num_nodes_ (num_nodes), 
    dimension_ (0)
{
    if (iterator_x != IT (0)) dimension_++;
    if (iterator_y != IT (0)) dimension_++;
    if (iterator_z != IT (0)) dimension_++;

    init_storage_ ();
    assign_data_ (iterator_x, iterator_y, iterator_z);
}

template <typename S>
void Coordinates<S>::swap (Coordinates<S>& rhs)
{
    std::swap (coordinates_, rhs.coordinates_);
    std::swap (num_nodes_, rhs.num_nodes_);
    std::swap (dimension_, rhs.dimension_);
}

template <typename S>
template <typename T>
void Coordinates<S>::take_data_from (std::vector<std::vector<T> >& rhs)
{
    const int dimensions = rhs.size ();
    for (int d = 0; d < dimensions; ++d)
        coordinates_ [d].assign (rhs [d].begin (), rhs [d].end ());

    dimension_ = coordinates_.size ();
    num_nodes_ = coordinates_ [0].size ();

    ASSERT (valid ());

}

template <typename S>
void Coordinates<S>::take_data_from (std::vector<std::vector<S> >& rhs)
{
    coordinates_.swap (rhs);

    dimension_ = coordinates_.size ();
    num_nodes_ = coordinates_ [0].size ();

    ASSERT (valid ());

}

template <typename S>
template <typename T>
Coordinates<S>* Coordinates<S>::build_from (std::vector<std::vector<T> >& rhs)
{
    Coordinates<S> *coords = new Coordinates<S>();
    
    coords->take_data_from (rhs);

    return coords;
}

template <typename S>
bool Coordinates<S>::valid () const
{

    bool result = (coordinates_.size () == dimension_);
    for (int d = 0; d < dimension_; ++d)
        result &= (coordinates_ [d].size () == num_nodes_);

    return result;

}


template <typename S>
void Coordinates<S>::to_stream (std::ostream& stream, bool verbose) const
{

    stream << "Coordinates:\n";
    stream << "  Nodes: " << num_nodes_ << "\n";
    stream << "  Dimensions: " << dimension_ << "\n";

    if (verbose)
    {
        for (int node = 0; node < num_nodes_ ; ++node)
        {
            stream << "Node " << node << ": " ;
            for (int dim = 0; dim < dimension_; ++dim)
            {
                stream << coordinates_ [dim] [node] << " ";
            }
            stream << "\n";
        }

    }
}


} // namespace

#endif
