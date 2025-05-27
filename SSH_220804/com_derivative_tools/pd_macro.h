#define PD_LOOP(i) for((i)=0;(i)<PD_VAR_END;(i)++)
#define PD_VAR0_LOOP(j) for((j)=0;(j)<PD_VAR0_END;(j)++)
#define PD_VAR1_LOOP(i) for((i)=0;(i)<PD_VAR1_END;(i)++)
#define PD_TOT_LOOP(j,i) PD_VAR0_LOOP(j)PD_VAR1_LOOP(i)
