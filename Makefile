SM=35
ARCH = -arch sm_$(SM)

#Build data
ifeq ($(OS),Windows_NT)
else
	GIT_DESCRIBE := "`git describe`"
	BUILD_DATE := "`date`"
	BUILD_SYSTEM := "`uname -a`"
	BUILD_PATH := "`pwd`"
endif
BUILD_DATA = -DGIT_DESCRIBE=\"$(GIT_DESCRIBE)\" -DBUILD_DATE=\"$(BUILD_DATE)\" -DBUILD_SYSTEM=\"$(BUILD_SYSTEM)\" -DBUILD_PATH=\"$(BUILD_PATH)\" -DBUILD_SM=\"$(SM)\"


source = heliosk.cu
headers = define.h ISO.h resample.h host.h voigt.h
all: heliosk hitran prepareExomol prepareExomolSuper
heliosk: $(source) $(headers)
	nvcc $(ARCH) --compiler-options -Wall -o heliosk $(source) $(BUILD_DATA)
#	nvcc -ccbin=g++-6 $(ARCH) -o heliosk $(source) $(BUILD_DATA)
#	nvcc --ptxas-options=-v $(ARCH) -o heliosk $(source) $(BUILD_DATA)
#	nvcc --compiler-options -Wall $(ARCH) -o heliosk $(source) $(BUILD_DATA)

hitran: hitran.cpp 
	g++ -Wall -Wextra -o hitran hitran.cpp
prepareExomol: prepareExomol.cpp ISO.h define.h
	g++ -Wall -Wextra -o prepareExomol prepareExomol.cpp
prepareExomolSuper: prepareExomolSuper.cpp ISO.h define.h
	g++ -Wall -Wextra -o prepareExomolSuper prepareExomolSuper.cpp

