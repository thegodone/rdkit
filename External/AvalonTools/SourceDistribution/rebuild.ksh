#!/usr/bin/ksh

export VERSION=` cat version.txt`

mkdir -p ../Binaries/classes
mkdir -p ../Binaries/archives
mkdir -p ../Binaries/$2/tomcat-lib
mkdir -p ../Binaries/$2/libraries

cp makefile ../Binaries/$2/makefile
echo cp $1.makefile ../Binaries/$2/platform.makefile

. ./$1_$2.defines

echo "###############################"
echo cat ./$1_$2.defines
cat ./$1_$2.defines
echo "###############################"

echo

echo "###############################"
pwd
ls -l
echo "###############################"

# Disable build of windows targets if not applicable
cp non_windows.makefile ../Binaries/$2/platform.makefile
if [[ "${WINDOWS_TARGETS}" = "yes" ]]; then
    cp windows.makefile ../Binaries/$2/platform.makefile
fi

cp license.txt ../Binaries/$2
cp readme_windows.txt ../Binaries/$2
cp readme_linux.txt ../Binaries/$2

cd java

echo
echo 'building javaclasses'
pwd
ls -l ../../StandardFiles/ThirdParty
make javaclasses
if [ $? -ne 0 ]; then
    echo 'building javaclasses failed'
    exit 1
fi
echo
echo 'building java stub header files'
make javastubs
if [ $? -ne 0 ]; then
    echo 'building Java stubs of C compilation failed'
    exit 1
fi
echo
echo 'building Tomcat WAR file'
ant
if [ $? -ne 0 ]; then
    echo 'building Tomcat WAR file failed'
    exit 1
fi

echo
echo 'building C artefacts'
cd -
cd ../Binaries/$2
make common
if [ $? -ne 0 ]; then
    echo 'building common targets failed'
    exit 1
fi
echo
echo 'building jni targets'
make jnitargets
if [ $? -ne 0 ]; then
    echo 'building jni target failed'
    exit 1
fi
echo
echo 'building Windows targets'
make wintargets
if [ $? -ne 0 ]; then
    echo 'building Windows targets failed'
    exit 1
fi

echo Copying files

pwd

cp ../archives/depictjni.jar tomcat-lib
cp ${LIB_PREFIX}JNIMatch.$LIB_EXT tomcat-lib
cp ${LIB_PREFIX}JNIDepict.$LIB_EXT tomcat-lib
cp avalon_jni_JNISmi2Mol.$LIB_EXT tomcat-lib
# make sure Linux systems also find the file
cp avalon_jni_JNISmi2Mol.$LIB_EXT tomcat-lib/${LIB_PREFIX}JNISmi2Mol.$LIB_EXT

cp ../archives/depictjni.jar libraries
cp ../archives/depict.jar libraries
cp ${LIB_PREFIX}JNIMatch.$LIB_EXT libraries
cp ${LIB_PREFIX}JNIDepict.$LIB_EXT libraries
cp avalon_jni_JNISmi2Mol.$LIB_EXT libraries
cp avalon_jni_JNISmi2Mol.$LIB_EXT libraries/${LIB_PREFIX}JNISmi2Mol.$LIB_EXT

cp ../archives/*.war .
cp ../../StandardFiles/checkfgs.* .
tar -cf ../AvalonToolkit_${VERSION}.$2.tar *
