#define NAN_ALLOC(a)  allocate(a); a = nan
#define SPVAL_ALLOC(a) allocate(a); a = spval
#define iSPVAL_ALLOC(a) allocate(a); a = ispval
#define ZERO_ALLOC(a) allocate(a); a = 0._r8
