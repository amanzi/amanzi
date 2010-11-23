#include <UnitTest++.h>
#include <stdexcept>

#include "../Coordinates.hh"

template <typename F>
void set (F v [], F a, F b, F c)
{
    v [0] = a;
    v [1] = b;
    v [2] = c;
}

template <typename F>
void set (std::vector<F> &v, F a, F b, F c)
{
    set (&v [0], a, b, c);
}

// Fixture base to hold a Coordinate object.
template <typename F>
struct Coords
{
    typedef Mesh_data::Coordinates<F> C;
    C* c;

    virtual ~Coords () { delete c; }
};

// Create an empty Coordinates<int> object.
struct No_Data : public Coords<int>
{
    No_Data ()  { c = new Mesh_data::Coordinates<int>(10,3);  }
};

// Create a coordinates object corresponding to the mesh in STK_mesh/test/Example_Mesh.hh
struct Data : public Coords<double>
{

    Data ()
    {
        std::vector<std::vector<double> > coord_data (3, std::vector<double>(20));
        for (int node = 0; node < 20; ++node)
        {
            const int x_plane  = node / 4;
            const int in_plane = node % 4;
            const int y_plane = int ((in_plane == 1) || (in_plane == 2));
            const int z_plane = int ((in_plane == 2) || (in_plane == 3));

            coord_data [0] [node] = double (x_plane);
            coord_data [1] [node] = double (y_plane);
            coord_data [2] [node] = double (z_plane);
        }
        c = Mesh_data::Coordinates<double>::build_from (coord_data);
    }

};


SUITE (Coordinates)
{

    TEST (Construction)
    {
        const int num_nodes = 10;

        std::vector<int> d (num_nodes*3,0);
        for (int i=0; i<d.size (); ++i)
            d [i] = i;

        Mesh_data::Coordinates<int> c1 (num_nodes, &d [0]);
        CHECK_EQUAL (c1.dimension (), 1);
        CHECK_EQUAL (c1.nodes (), num_nodes);
        CHECK_ARRAY_EQUAL (c1.get_coordinate_pointer (0), &d[0], 10);

        Mesh_data::Coordinates<int> c2 (num_nodes, &d [0], &d [10]);
        CHECK_EQUAL (c2.dimension (), 2);
        CHECK_EQUAL (c1.nodes (), num_nodes);
        CHECK_ARRAY_EQUAL (c2.get_coordinate_pointer (0), &d [0 ], 10);
        CHECK_ARRAY_EQUAL (c2.get_coordinate_pointer (1), &d [10], 10);

        Mesh_data::Coordinates<int> c3 (num_nodes, &d [0], &d [10], &d [20]);
        CHECK_EQUAL (c3.dimension (), 3);
        CHECK_EQUAL (c1.nodes (), num_nodes);
        CHECK_ARRAY_EQUAL (c3.get_coordinate_pointer (0), &d [0 ], 10);
        CHECK_ARRAY_EQUAL (c3.get_coordinate_pointer (1), &d [10], 10);
        CHECK_ARRAY_EQUAL (c3.get_coordinate_pointer (2), &d [20], 10);

    }

    TEST (Size)
    {
        
        Mesh_data::Coordinates<double> d(10, 2);

        CHECK_EQUAL (d.nodes (), 10);
        CHECK_EQUAL (d.dimension (), 2);

        Mesh_data::Coordinates<float> f(100,1);
        
        CHECK_EQUAL (f.nodes (), 100);
        CHECK_EQUAL (f.dimension (), 1);

        Mesh_data::Coordinates<int> i(1000, 3);
        
        CHECK_EQUAL (i.nodes (), 1000);
        CHECK_EQUAL (i.dimension (), 3);

    }

    TEST_FIXTURE (No_Data, Storage)
    {
        (*c) (2,0) = 10;
        CHECK_EQUAL ((*c) (2,0), 10);
    }

    TEST_FIXTURE (No_Data, Pointer_assignment)
    {

        C::pointer_type p = c->get_coordinate_pointer (0);

        p [5] = 10;
        CHECK_EQUAL ((*c) (5,0), 10);
        
        p = c->get_coordinate_pointer (2);
    }

    TEST_FIXTURE (No_Data, Bounds_checking)
    {
        CHECK_THROW ((*c) (10, 1), DBC::Assertion);
        CHECK_THROW ((*c) (-1, 1), DBC::Assertion);
        CHECK_THROW ((*c) (5,  3), DBC::Assertion);
        CHECK_THROW ((*c) (5, -1), DBC::Assertion);
    }


    TEST_FIXTURE (No_Data, Assignment)
    {
        std::vector<std::vector<int> > d(3, std::vector<int>(10));
        (*c).take_data_from (d);
    }

    TEST_FIXTURE (Data, Retrevial)
    {
        double coords[3];
        double values[3];

        set (values, 0.0, 0.0, 0.0);
        (*c) (0, coords);
        CHECK_ARRAY_EQUAL (coords, values, 3);

        set (values, 4.0, 1.0, 1.0);
        (*c) (18, coords);
        CHECK_ARRAY_EQUAL (coords, values, 3);
    }
    

}
