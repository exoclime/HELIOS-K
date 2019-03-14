SM=35
ARCH = -arch sm_$(SM)


GIT_DESCRIBE := "`git describe`"
BUILD_DATE := "`date`"
BUILD_SYSTEM := "`uname -a`"
BUILD_PATH := "`pwd`"
BUILD_DATA = -DGIT_DESCRIBE=\"$(GIT_DESCRIBE)\" -DBUILD_DATE=\"$(BUILD_DATE)\" -DBUILD_SYSTEM=\"$(BUILD_SYSTEM)\" -DBUILD_PATH=\"$(BUILD_PATH)\" -DBUILD_SM=\"$(SM)\"


source = heliosk.cu
headers = define.h ISO.h resample.h host.h voigt.h
all: heliosk prepare prepareExomol
heliosk: $(source) $(headers)
	nvcc $(ARCH) --compiler-options -Wall -o heliosk $(source) $(BUILD_DATA)
#	nvcc -ccbin=g++-6 $(ARCH) -o heliosk $(source) $(BUILD_DATA)
#	nvcc --ptxas-options=-v $(ARCH) -o heliosk $(source) $(BUILD_DATA)
#	nvcc --compiler-options -Wall $(ARCH) -o heliosk $(source) $(BUILD_DATA)

prepare: prepare.cpp ISO.h define.h
	g++ -o prepare prepare.cpp
prepareExomol: prepareExomol.cpp ISO.h define.h
	g++ -o prepareExomol prepareExomol.cpp

