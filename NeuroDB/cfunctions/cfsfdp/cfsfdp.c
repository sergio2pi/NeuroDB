//compiling with : gcc -fPIC -o cfsfdp.so -shared cfsfdp.c -I/usr/include/pgsql -L/usr/lib/postgresql92/lib64 -lpq

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <libpq-fe.h>
#include <string.h>
#include <math.h>


int compare(const void *_a, const void *_b) {
    int *a, *b;

    a = (int *) _a;
    b = (int *) _b;

    return (*a - *b);
}

float get_distance(float x1, float y1, float z1, float x2, float y2, float z2)
{
    double tmp = pow(x1 - x2, 2) + pow(y1 - y2, 2) + pow(z1 - z2, 2);
    return pow(tmp, 0.5);
}


float get_distance_from_res(PGresult *res, int i, int j)
{
    float tmp;
    float x1, x2, y1, y2, z1, z2;
    sscanf(PQgetvalue(res, i, 0),"%f",&x1);
    sscanf(PQgetvalue(res, i, 1),"%f",&y1);
    sscanf(PQgetvalue(res, i, 2),"%f",&z1);
    sscanf(PQgetvalue(res, j, 0),"%f",&x2);
    sscanf(PQgetvalue(res, j, 1),"%f",&y2);
    sscanf(PQgetvalue(res, j, 2),"%f",&z2);
    tmp = get_distance(x1, y1, z1, x2, y2, z2);

    return tmp;

}

float get_dc(char connect[150], char id_block[20], float percent)
{
    char query[200];
    strcpy(query,"SELECT spike.p1, spike.p2, spike.p3, spike.id from SPIKE JOIN  segment ON id_segment = segment.id WHERE segment.id_block = ");
    strcat(query, id_block);

    PGconn          *conn;
    PGresult        *res;
    int             rec_count;

    int n = 0, i, j;
    float dc = 0.0;
    float* distances = NULL;
    int position;

    conn = PQconnectdb(connect);

    if (PQstatus(conn) == CONNECTION_BAD) {
        puts("We were unable to connect to the database");
        return 0.0;
    }
    res = PQexec(conn,query);

    if (PQresultStatus(res) != PGRES_TUPLES_OK) {
        puts("We did not get any data!");
        return 0.0;
    }

    rec_count = PQntuples(res);
    distances = (float*)malloc(sizeof(float)*rec_count*rec_count);

    for (i = 0, n = 0; i < rec_count; i++) {
          for (j = 0; j < rec_count; j++, n++) {
                distances[n] = get_distance_from_res(res, i, j);
          }
      }

    //There are rec_count zeros, because distances from itselfs
    position = rec_count*percent/100 -1;//position = rec_count + 2*rec_count*percent/100 -1;
    qsort(distances, rec_count, sizeof(float), &compare);

    dc = distances[position];

    PQclear(res);
    PQfinish(conn);
    return dc;
}

int get_local_density(char connect[150], char id_block[20], float dc, double* local_density, char kernel[20])
//TODO: Parameter size
{

    PGconn          *conn;
    PGresult        *res;
    int             rec_count;

    char query[200];
    int i, j;
    double distance;

    if (strcmp(kernel, "gaussian") && strcmp(kernel, "cutoff"))
    {
        puts("Kernel parameter must be \"gaussian\" or \"cutoff\".");
        return 1;
    }

    strcpy(query,"SELECT spike.p1, spike.p2, spike.p3 from SPIKE JOIN  segment ON id_segment = segment.id WHERE segment.id_block = ");
    strcat(query, id_block);

    conn = PQconnectdb(connect);

    if (PQstatus(conn) == CONNECTION_BAD) {
        puts("We were unable to connect to the database");
        return 1;
    }


    res = PQexec(conn,query);

    if (PQresultStatus(res) != PGRES_TUPLES_OK) {
        puts("We did not get any data!");
        return 2;
    }

    rec_count = PQntuples(res);

    for (i = 0; i < rec_count; i++) local_density[i] = 0;

    if (!strcmp(kernel, "cutoff") )
    {
        for (i = 0; i < rec_count; i++) {
              for (j = i+1; j < rec_count; j++) {
                    distance = get_distance_from_res(res, i, j);
                    if ( (distance - dc) < 0 )
                    {
                        local_density[i] += 1;
                        local_density[j] += 1;
                    }
              }
          }
    }

    if (!strcmp(kernel, "gaussian") )
    {
        for (i = 0; i < rec_count; i++) {
              for (j = i+1; j < rec_count; j++) {
                    distance = get_distance_from_res(res, i, j);
                    local_density[i] = local_density[i] + exp(-(distance/dc)*(distance/dc));
                    local_density[j] = local_density[j] + exp(-(distance/dc)*(distance/dc));
              }
          }
    }

    PQclear(res);
    PQfinish(conn);

    return 0;
}


int get_distance_to_higher_density(char connect[], char id_block[], double* rho, double* delta, int size){

    PGconn          *conn;
    PGresult        *res;
    int             rec_count;

    int nSamples = size;
    double dist;
    double tmp;
    char query[200];
    int i, j, k, flag;

    strcpy(query,"SELECT spike.p1, spike.p2, spike.p3 from SPIKE JOIN  segment ON id_segment = segment.id WHERE segment.id_block = ");
    strcat(query, id_block);

    conn = PQconnectdb(connect);

    if (PQstatus(conn) == CONNECTION_BAD) {
        puts("We were unable to connect to the database");
        return 1;
    }

    res = PQexec(conn,query);

    if (PQresultStatus(res) != PGRES_TUPLES_OK) {
        puts("We did not get any data!");
        return 2;
    }

    rec_count = PQntuples(res);

    for (i = 0; i < rec_count; i++)
        delta[i] = 0;

    for(i = 0; i < nSamples; i++){
        dist = 0.0;
        flag = 0;
        for(j = 0; j < nSamples; j++){
            if(i == j) continue;
            if(rho[j] > rho[i]){

                tmp = get_distance_from_res(res, i, j);

                if(!flag){
                    dist = tmp;
                    flag = 1;
                }else dist = tmp < dist ? tmp : dist;
            }
        }
        if(!flag){
            for(k = 0; k < nSamples; k++){
                tmp = get_distance_from_res(res, i, k);
                dist = tmp > dist ? tmp : dist;
            }
        }
        delta[i] = dist;
    }

    PQclear(res);
    PQfinish(conn);

    return 0;
}

double mean(double v[], int n)
{
    int i;
    float m = 0;

    for(i=0; i<n; i++) m += v[i];

    return m/n;
}

double* get_centers_cluster_dp(char connect[], char id_block[])
{
    PGconn          *conn;
    PGresult        *res;
    int             rec_count;

    char query[200];
    float dc;
    double* local_density;
    double* distance_to_higher_density;
    double* gamma;
    float meanf;
    double* centers = NULL;
    int k=0, i;

    conn = PQconnectdb(connect);

    strcpy(query,"SELECT spike.id, spike.p1, spike.p2, spike.p3 from SPIKE JOIN  segment ON id_segment = segment.id WHERE segment.id_block = ");
    strcat(query, id_block);

    if (PQstatus(conn) == CONNECTION_BAD) {
        puts("cluster_dp: We were unable to connect to the database");
        return NULL;
    }

    res = PQexec(conn,query);
    rec_count = PQntuples(res);
    local_density = (double*)calloc(rec_count,sizeof(double));
    distance_to_higher_density = (double*)calloc(rec_count,sizeof(double));
    gamma = (double*)calloc(rec_count,sizeof(double));

    dc = get_dc("dbname=demo host=192.168.2.2 user=postgres password=postgres", "54", 2);

    get_local_density("dbname=demo host=192.168.2.2 user=postgres password=postgres", "54", dc, local_density, "cutoff");
    get_distance_to_higher_density("dbname=demo host=192.168.2.2 user=postgres password=postgres", "54",local_density, distance_to_higher_density, 1026);

    // gamma is rho*delta
    for(i=0; i<rec_count; i++)
    {
        gamma[i] = local_density[i]*distance_to_higher_density[i];
    }

    //points > 5*mean are centers of cluster
    meanf = mean(gamma, rec_count);
    centers = (double*)malloc(sizeof(double));
    centers[0] = 0;
    for(i=0; i<rec_count; i++)
    {
        if (gamma[i] > 5*meanf)
        {
            centers = (double*)realloc(centers, sizeof(double)*(centers[0]+1));
            sscanf(PQgetvalue(res, i, 0),"%lf",&centers[(int)centers[0]+1]);
            printf("%lf\n",centers[(int)centers[0]+1]);
            centers[0]++;
        }
    }

    // centers contains id of spike centers
    return centers;
}

double** get_cluster_dp(char connect[], char id_block[], double* centers)
{
    PGconn          *conn;
    PGresult        *res;
    int             rec_count;

    char query[200];

    conn = PQconnectdb(connect);

    strcpy(query,"SELECT spike.id, spike.p1, spike.p2, spike.p3 from SPIKE JOIN  segment ON id_segment = segment.id WHERE segment.id_block = ");
    strcat(query, id_block);

    if (PQstatus(conn) == CONNECTION_BAD) {
        puts("cluster_dp: We were unable to connect to the database");
        return NULL;
    }

    res = PQexec(conn,query);
    rec_count = PQntuples(res);


}

int main(int argc, char* argv[])
{
//float dc = 0;
    int i;
    double* local_density = (double*)calloc(1026,sizeof(double));
    double* distance_to_higher_density = (double*)calloc(1026,sizeof(double));
    //double* gamma = (double*)calloc(1026,sizeof(double));
    double* centers = NULL;

    float dc = get_dc("dbname=demo host=192.168.2.2 user=postgres password=postgres", "54", 2);

    get_local_density("dbname=demo host=192.168.2.2 user=postgres password=postgres", "54", dc, local_density, "gaussian");

    get_distance_to_higher_density("dbname=demo host=192.168.2.2 user=postgres password=postgres", "54",local_density, distance_to_higher_density, 1026);

    centers = get_centers_cluster_dp("dbname=demo host=192.168.2.2 user=postgres password=postgres", "54");

    return 1;
}
