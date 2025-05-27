#define ARRAY_LOOP(i) for(i = 0; i < g_mtx_size; i++)
#define MTX_LOOP(j,i) ARRAY_LOOP(j)ARRAY_LOOP(i)

#define NX1_LOOP(i) for(i = 0; i < g_mtx_nx1; i++)
#define NX2_LOOP(i) for(i = 0; i < g_mtx_nx2; i++)
#define NX21_LOOP(j,i) NX2_LOOP(j)NX1_LOOP(i)

