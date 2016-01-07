SM=35
ARCH = -arch sm_$(SM)
source = heliosk.cu
headers = define.h ISO.h resample.h host.h voigt.h
heliosk: $(source) $(headers)
	nvcc $(ARCH) -o heliosk $(source)
#	nvcc --ptxas-options=-v $(ARCH) -o heliosk $(source)
#	nvcc --compiler-options -Wall $(ARCH) -o heliosk $(source)
