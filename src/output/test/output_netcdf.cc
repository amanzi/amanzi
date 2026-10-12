#include <mpi.h>
#include <netcdf.h>
#include <netcdf_par.h>

#include <stdio.h>
#include <stdlib.h>

#define NC_CHECK(call)                                      \
    do {                                                    \
        int status = (call);                                \
        if (status != NC_NOERR) {                           \
            fprintf(stderr, "[%d] %s:%d: %s\n",             \
                    rank, __FILE__, __LINE__,               \
                    nc_strerror(status));                  \
            MPI_Abort(MPI_COMM_WORLD, status);              \
        }                                                   \
    } while (0)


int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    int rank, nranks;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nranks);

    /* -----------------------------------------------------
     * Global grid
     * ----------------------------------------------------- */

    const size_t NLAT = 4000;
    const size_t NLON = 8000;

    /*
     * Simple block decomposition in latitude.
     *
     * In a real code these would presumably come from
     * the model's existing domain decomposition.
     */
    size_t base = NLAT / nranks;
    size_t rem  = NLAT % nranks;

    size_t lat_count =
        base + ((size_t)rank < rem ? 1 : 0);

    size_t lat_start =
        rank * base +
        ((size_t)rank < rem ? (size_t)rank : rem);


    /* Local temperature field for one timestep. */
    double *temperature =
      (double*)malloc(lat_count * NLON * sizeof(double));

    if (!temperature) {
        fprintf(stderr, "allocation failed\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }


    /* -----------------------------------------------------
     * Create the NetCDF file collectively
     * ----------------------------------------------------- */

    int ncid;

    NC_CHECK(
        nc_create_par("output.nc",
                      NC_NETCDF4 | NC_CLOBBER,
                      MPI_COMM_WORLD,
                      MPI_INFO_NULL,
                      &ncid)
    );


    /* -----------------------------------------------------
     * Define dimensions -- COLLECTIVE
     * ----------------------------------------------------- */

    int time_dim, lat_dim, lon_dim;

    NC_CHECK(nc_def_dim(ncid, "time", NC_UNLIMITED, &time_dim));
    NC_CHECK(nc_def_dim(ncid, "lat",  NLAT,         &lat_dim));
    NC_CHECK(nc_def_dim(ncid, "lon",  NLON,         &lon_dim));


    /* -----------------------------------------------------
     * Define variables -- COLLECTIVE
     * ----------------------------------------------------- */

    int time_var;
    int lat_var;
    int lon_var;
    int temperature_var;
    int total_water_var;

    int time_dims[1] = { time_dim };
    int lat_dims[1]  = { lat_dim };
    int lon_dims[1]  = { lon_dim };

    int temperature_dims[3] = {
        time_dim, lat_dim, lon_dim
    };

    int total_water_dims[1] = { time_dim };

    NC_CHECK(nc_def_var(ncid, "time",
                        NC_DOUBLE, 1, time_dims, &time_var));

    NC_CHECK(nc_def_var(ncid, "lat",
                        NC_DOUBLE, 1, lat_dims, &lat_var));

    NC_CHECK(nc_def_var(ncid, "lon",
                        NC_DOUBLE, 1, lon_dims, &lon_var));

    NC_CHECK(nc_def_var(ncid, "temperature",
                        NC_DOUBLE, 3, temperature_dims,
                        &temperature_var));

    // NC_CHECK(nc_def_var(ncid, "total_water",
    //                     NC_DOUBLE, 1, total_water_dims,
    //                     &total_water_var));


    /* -----------------------------------------------------
     * Attributes -- also collective
     * ----------------------------------------------------- */

    NC_CHECK(nc_put_att_text(
        ncid, time_var, "units",
        20, "seconds since 2000-01-01"));

    NC_CHECK(nc_put_att_text(
        ncid, lat_var, "units",
        7, "degrees"));

    NC_CHECK(nc_put_att_text(
        ncid, lon_var, "units",
        7, "degrees"));

    NC_CHECK(nc_put_att_text(
        ncid, temperature_var, "units",
        1, "K"));

    // NC_CHECK(nc_put_att_text(
    //     ncid, total_water_var, "units",
    //     1, "m"));


    /*
     * Temperature is the large distributed variable.
     *
     * Make ONLY this variable collective.
     */
    NC_CHECK(nc_var_par_access(
        ncid, temperature_var, NC_COLLECTIVE));


    /*
     * The other variables remain independent.
     *
     * NC_INDEPENDENT is the default, so this is not strictly
     * necessary, but making it explicit makes the intent clear.
     */
    NC_CHECK(nc_var_par_access(
        ncid, time_var, NC_INDEPENDENT));

    NC_CHECK(nc_var_par_access(
        ncid, lat_var, NC_INDEPENDENT));

    NC_CHECK(nc_var_par_access(
        ncid, lon_var, NC_INDEPENDENT));

    // NC_CHECK(nc_var_par_access(
    //     ncid, total_water_var, NC_INDEPENDENT));


    /* -----------------------------------------------------
     * Finish definition mode -- COLLECTIVE
     * ----------------------------------------------------- */

    NC_CHECK(nc_enddef(ncid));


    /* -----------------------------------------------------
     * Write static coordinates.
     *
     * Only rank 0 needs to do this.
     * ----------------------------------------------------- */

    if (rank == 0) {

      double *lat = (double*)malloc(NLAT * sizeof(double));
      double *lon = (double*)malloc(NLON * sizeof(double));

        for (size_t i = 0; i < NLAT; ++i)
            lat[i] = -90.0 + 180.0 * i / (NLAT - 1);

        for (size_t j = 0; j < NLON; ++j)
            lon[j] = -180.0 + 360.0 * j / NLON;

        NC_CHECK(nc_put_var_double(
            ncid, lat_var, lat));

        NC_CHECK(nc_put_var_double(
            ncid, lon_var, lon));

        free(lat);
        free(lon);
    }

    /*
     * Not strictly required for correctness of the NetCDF
     * variable operations below, but this makes the phase
     * transition explicit.
     */
    MPI_Barrier(MPI_COMM_WORLD);


    /* -----------------------------------------------------
     * Time loop
     * ----------------------------------------------------- */

    const int NT = 1000;

    for (int t = 0; t < NT; ++t) {

        /*
         * -----------------------------------------------
         * Compute local temperature[lat_count][NLON]
         * -----------------------------------------------
         */

        for (size_t i = 0; i < lat_count; ++i) {
            for (size_t j = 0; j < NLON; ++j) {

                /*
                 * Placeholder for model data.
                 */
                temperature[i * NLON + j] =
                    273.15 + rank + 0.001 * t;
            }
        }


        /*
         * -----------------------------------------------
         * Compute total water.
         *
         * Each rank has its local contribution.
         * -----------------------------------------------
         */

        double local_water = 0.0;

        for (size_t i = 0; i < lat_count; ++i) {
            for (size_t j = 0; j < NLON; ++j) {
                local_water += temperature[i * NLON + j];
            }
        }

        double total_water;

        MPI_Allreduce(&local_water,
                      &total_water,
                      1,
                      MPI_DOUBLE,
                      MPI_SUM,
                      MPI_COMM_WORLD);


        /*
         * -----------------------------------------------
         * 1. Temperature: COLLECTIVE
         * -----------------------------------------------
         *
         * Every rank participates.
         *
         * Each rank writes:
         *
         *   time = t
         *   lat  = [lat_start, lat_start + lat_count)
         *   lon  = [0, NLON)
         */

        size_t temp_start[3] = {
            (size_t)t,
            lat_start,
            0
        };

        size_t temp_count[3] = {
            1,
            lat_count,
            NLON
        };

        NC_CHECK(nc_put_vara_double(
            ncid,
            temperature_var,
            temp_start,
            temp_count,
            temperature));


        /*
         * -----------------------------------------------
         * 2. total_water: INDEPENDENT
         * -----------------------------------------------
         *
         * total_water is replicated because of Allreduce().
         *
         * There is no reason for all ranks to write it.
         *
         * Rank 0 writes the single element.
         */

        // if (rank == 0) {

        //     size_t scalar_start[1] = { (size_t)t };
        //     size_t scalar_count[1] = { 1 };

        //     NC_CHECK(nc_put_vara_double(
        //         ncid,
        //         total_water_var,
        //         scalar_start,
        //         scalar_count,
        //         &total_water));
        // }


        /*
         * -----------------------------------------------
         * 3. time: INDEPENDENT
         * -----------------------------------------------
         *
         * Again, rank 0 owns this tiny write.
         */

        if (rank == 0) {

            size_t time_start[1] = { (size_t)t };
            size_t time_count[1] = { 1 };

            double current_time = 3600.0 * t;

            NC_CHECK(nc_put_vara_double(
                ncid,
                time_var,
                time_start,
                time_count,
                &current_time));
        }
    }


    /* -----------------------------------------------------
     * Close -- collective
     * ----------------------------------------------------- */

    NC_CHECK(nc_close(ncid));

    free(temperature);

    MPI_Finalize();

    return 0;
}
