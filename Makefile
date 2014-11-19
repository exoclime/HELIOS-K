SM=35
ARCH = -arch sm_$(SM)
source = kcalc.cu
headers = define.h ISO.h
kcalc: $(source) $(headers)
	nvcc $(ARCH) -o kcalc $(source)
