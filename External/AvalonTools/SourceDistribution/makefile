# define variables if not predefined
LIB_EXT ?= so
OBJ_EXT ?= o
CC ?= cc
LD ?= cc -lm
BITS ?= 32
## 
include platform.makefile
#
# standard targets
common : matchtest${EXE} canonizer${EXE} smi2mol${EXE} mol2smi${EXE} mol2tbl${EXE} struchk${EXE}
#
# Java class files, use jni directory as a proxy to trigger compilation of all Java sources
# JAVA_LIBS = lib/log4j-1.2.8.jar;lib/commons-logging-1.1.jar;lib/servlet.jar
JAVA_LIBS = ../../Binaries/lib/log4j-1.2.8.jar
#
javaclasses :
	javac -cp "../../Binaries/classes;${JAVA_LIBS}" -d ../../Binaries/classes -sourcepath . jni/JNIDepict.java
	javac -Xlint:unchecked -cp "../../Binaries/classes;${JAVA_LIBS}" -d ../../Binaries/classes -sourcepath . jni/JNIMatch.java
	javac -cp "../../Binaries/classes;${JAVA_LIBS}" -d ../../Binaries/classes -sourcepath . jni/JNIWinTools.java
	javac -cp "../../Binaries/classes;${JAVA_LIBS}" -d ../../Binaries/classes -sourcepath . jni/UnixTools.java
	javac -cp "../../Binaries/classes;${JAVA_LIBS}" -d ../../Binaries/classes -sourcepath . avalon/jni/JNISmi2Mol.java
	sed "s/#VERSION_NUMBER#/${VERSION}/" <avalon/jsp/servlets/Mol2ImageServlet.java >avalon/jsp/servlets/Mol2ImageServlet.java_versioned
	cp avalon/jsp/servlets/Mol2ImageServlet.java_versioned avalon/jsp/servlets/Mol2ImageServlet.java
	javac -Xlint:deprecation,unchecked -cp "../../Binaries/classes;${JAVA_LIBS}" -d ../../Binaries/classes -sourcepath . avalon/jsp/servlets/Mol2ImageServlet.java
# Create header files for C code
javastubs :
	javah -classpath "../../Binaries/classes" -d ../jni jni.JNIDepict
	javah -classpath "../../Binaries/classes" -d ../jni jni.JNIMatch
	javah -classpath "../../Binaries/classes" -d ../jni jni.JNIWinTools
	javah -classpath "../../Binaries/classes" -d ../jni avalon.jni.JNISmi2Mol
#
# JNI targets
jnitargets : ${LIB_PREFIX}JNIDepict.${LIB_EXT} avalon_jni_JNISmi2Mol.${LIB_EXT} ${LIB_PREFIX}JNIMatch.${LIB_EXT} librohde_clib.${LIB_EXT}
#
#
JDK_INCLUDE ?= c:/jdk6/include
JDK_MACHINE_INCLUDE ?= c:/jdk6/include/win32

JNIDIR = ../../SourceDistribution/jni
WINDIR = ../../SourceDistribution/mswindows
COMMONDIR = ../../SourceDistribution/common
#
# This is a useful URL describing among other calling convention and name mangling issues the kill-at option
# http://www.willus.com/mingw/yongweiwu_stdcall.html
toolkit${BITS}.${LIB_EXT} : \
	       canonizer.${OBJ_EXT} \
	       casutils.${OBJ_EXT} \
	       denormal.${OBJ_EXT} \
	       forio.${OBJ_EXT} \
	       geometry.${OBJ_EXT} \
	       graph.${OBJ_EXT} \
	       hashcode.${OBJ_EXT} \
	       layout.${OBJ_EXT} \
	       local.${OBJ_EXT} \
	       perceive.${OBJ_EXT} \
	       rtutils.${OBJ_EXT} \
	       reaccsio.${OBJ_EXT} \
	       set.${OBJ_EXT} \
	       smi2mol.${OBJ_EXT} pattern.${OBJ_EXT} \
	       shortcut.${OBJ_EXT} \
	       ssmatch.${OBJ_EXT} \
	       stereo.${OBJ_EXT} \
	       symbol_lists.${OBJ_EXT} \
	       symboltable.${OBJ_EXT} \
	       utilities.${OBJ_EXT} \
	       toolkit.${OBJ_EXT}
	${LD} -shared -o toolkit${BITS}.${LIB_EXT} \
	--verbose \
    ${KILL_AT} \
    --subsystem,windows \
    -mwindows \
    canonizer.${OBJ_EXT} \
    casutils.${OBJ_EXT} \
    denormal.${OBJ_EXT} \
    forio.${OBJ_EXT} \
    geometry.${OBJ_EXT} \
    graph.${OBJ_EXT} \
    hashcode.${OBJ_EXT} \
    layout.${OBJ_EXT} \
    local.${OBJ_EXT} \
    perceive.${OBJ_EXT} \
    reaccsio.${OBJ_EXT} \
    rtutils.${OBJ_EXT} \
    set.${OBJ_EXT} \
    smi2mol.${OBJ_EXT} pattern.${OBJ_EXT} \
    shortcut.${OBJ_EXT} \
    ssmatch.${OBJ_EXT} \
    stereo.${OBJ_EXT} \
    symbol_lists.${OBJ_EXT} \
    symboltable.${OBJ_EXT} \
    utilities.${OBJ_EXT} \
    toolkit.${OBJ_EXT} \
	${LD_OPTS}
#
depict${BITS}.${LIB_EXT} : casutils.${OBJ_EXT} denormal.${OBJ_EXT} forio.${OBJ_EXT} geometry.${OBJ_EXT} graph.${OBJ_EXT} \
	       layout.${OBJ_EXT} local.${OBJ_EXT} perceive.${OBJ_EXT} reaccsio.${OBJ_EXT} rtutils.${OBJ_EXT} \
	       set.${OBJ_EXT} sketch.${OBJ_EXT} smi2mol.${OBJ_EXT} pattern.${OBJ_EXT} shortcut.${OBJ_EXT} symboltable.${OBJ_EXT} stereo.${OBJ_EXT} symbol_lists.${OBJ_EXT} \
               utilities.${OBJ_EXT} hashcode.${OBJ_EXT} ssmatch.${OBJ_EXT}  \
	       depict.${OBJ_EXT}
	${CC} -shared -o depict${BITS}.${LIB_EXT} \
    ${KILL_AT} \
    -Wl,--subsystem,windows \
    -mwindows \
    casutils.${OBJ_EXT} denormal.${OBJ_EXT} forio.${OBJ_EXT} geometry.${OBJ_EXT} graph.${OBJ_EXT} \
    layout.${OBJ_EXT} local.${OBJ_EXT} perceive.${OBJ_EXT} reaccsio.${OBJ_EXT} rtutils.${OBJ_EXT} \
    set.${OBJ_EXT} sketch.${OBJ_EXT} smi2mol.${OBJ_EXT} pattern.${OBJ_EXT} shortcut.${OBJ_EXT} symboltable.${OBJ_EXT} stereo.${OBJ_EXT} symbol_lists.${OBJ_EXT} \
    hashcode.${OBJ_EXT} ssmatch.${OBJ_EXT} \
    utilities.${OBJ_EXT} \
    depict.${OBJ_EXT}
#
canonizer${EXE} : canonizer_main.${OBJ_EXT} \
    casutils.${OBJ_EXT} denormal.${OBJ_EXT} forio.${OBJ_EXT} geometry.${OBJ_EXT} \
    graph.${OBJ_EXT} layout.${OBJ_EXT} perceive.${OBJ_EXT} \
    reaccsio.${OBJ_EXT} rtutils.${OBJ_EXT} set.${OBJ_EXT} smi2mol.${OBJ_EXT} pattern.${OBJ_EXT} shortcut.${OBJ_EXT} \
    stereo.${OBJ_EXT} symbol_lists.${OBJ_EXT} symboltable.${OBJ_EXT} utilities.${OBJ_EXT} \
    local.${OBJ_EXT} hashcode.${OBJ_EXT} ssmatch.${OBJ_EXT}
	${LD} -o canonizer${EXE} canonizer_main.${OBJ_EXT} \
	    casutils.${OBJ_EXT} denormal.${OBJ_EXT} forio.${OBJ_EXT} geometry.${OBJ_EXT} \
	    graph.${OBJ_EXT} layout.${OBJ_EXT} perceive.${OBJ_EXT} \
	    reaccsio.${OBJ_EXT} rtutils.${OBJ_EXT} set.${OBJ_EXT} smi2mol.${OBJ_EXT} pattern.${OBJ_EXT} shortcut.${OBJ_EXT} \
	    stereo.${OBJ_EXT} symbol_lists.${OBJ_EXT} symboltable.${OBJ_EXT} utilities.${OBJ_EXT} \
	    local.${OBJ_EXT} hashcode.${OBJ_EXT} ssmatch.${OBJ_EXT} \
	    ${LD_OPTS}
#
smi2mol${EXE} : smi2mol_main.${OBJ_EXT} \
    casutils.${OBJ_EXT} denormal.${OBJ_EXT} didepict.${OBJ_EXT} forio.${OBJ_EXT} geometry.${OBJ_EXT} \
    graph.${OBJ_EXT} layout.${OBJ_EXT} perceive.${OBJ_EXT} \
    reaccsio.${OBJ_EXT} rtutils.${OBJ_EXT} set.${OBJ_EXT} shortcut.${OBJ_EXT} \
    stereo.${OBJ_EXT} symbol_lists.${OBJ_EXT} symboltable.${OBJ_EXT} utilities.${OBJ_EXT} \
    local.${OBJ_EXT} hashcode.${OBJ_EXT} ssmatch.${OBJ_EXT} 
	${LD} -o smi2mol${EXE} smi2mol_main.${OBJ_EXT} pattern.${OBJ_EXT} shortcut.${OBJ_EXT} \
	    casutils.${OBJ_EXT} denormal.${OBJ_EXT} didepict.${OBJ_EXT} forio.${OBJ_EXT} geometry.${OBJ_EXT} \
	    graph.${OBJ_EXT} layout.${OBJ_EXT} perceive.${OBJ_EXT} \
	    reaccsio.${OBJ_EXT} rtutils.${OBJ_EXT} set.${OBJ_EXT} \
	    stereo.${OBJ_EXT} symbol_lists.${OBJ_EXT} symboltable.${OBJ_EXT} utilities.${OBJ_EXT} \
	    local.${OBJ_EXT} hashcode.${OBJ_EXT} ssmatch.${OBJ_EXT} \
	    ${LD_OPTS}
#
mol2smi${EXE} : smi2mol${EXE}
	cp smi2mol${EXE} mol2smi${EXE}
#
mol2tbl${EXE} : smi2mol${EXE}
	cp smi2mol${EXE} mol2tbl${EXE}
#
matchtest${EXE} : matchtest.${OBJ_EXT} \
    casutils.${OBJ_EXT} \
    denormal.${OBJ_EXT} \
    forio.${OBJ_EXT} \
    geometry.${OBJ_EXT} \
    graph.${OBJ_EXT} \
    layout.${OBJ_EXT} \
    perceive.${OBJ_EXT} \
    reaccsio.${OBJ_EXT} \
    rtutils.${OBJ_EXT} \
    set.${OBJ_EXT} \
    smi2mol.${OBJ_EXT} \
    pattern.${OBJ_EXT} \
    shortcut.${OBJ_EXT} \
    stereo.${OBJ_EXT} \
    symbol_lists.${OBJ_EXT} \
    symboltable.${OBJ_EXT} \
    utilities.${OBJ_EXT} \
    local.${OBJ_EXT} \
    hashcode.${OBJ_EXT} \
    ssmatch.${OBJ_EXT}
	${LD} -o matchtest${EXE} matchtest.${OBJ_EXT} \
	    casutils.${OBJ_EXT} \
	    denormal.${OBJ_EXT} \
	    forio.${OBJ_EXT} \
	    geometry.${OBJ_EXT} \
	    graph.${OBJ_EXT} \
	    layout.${OBJ_EXT} \
	    perceive.${OBJ_EXT} \
	    reaccsio.${OBJ_EXT} \
	    rtutils.${OBJ_EXT} \
	    set.${OBJ_EXT} \
	    smi2mol.${OBJ_EXT} \
	    pattern.${OBJ_EXT} \
	    shortcut.${OBJ_EXT} \
	    stereo.${OBJ_EXT} \
	    symbol_lists.${OBJ_EXT} \
	    symboltable.${OBJ_EXT} \
	    utilities.${OBJ_EXT} \
	    local.${OBJ_EXT} \
	    hashcode.${OBJ_EXT} \
	    ssmatch.${OBJ_EXT} \
	    ${LD_OPTS}
	cp matchtest${EXE} mol2fp${EXE}
	cp matchtest${EXE} smi2fp${EXE}
# Special compiles for main() modules
#
canonizer_main.${OBJ_EXT} : ${COMMONDIR}/canonizer.c  ${COMMONDIR}/canonizer.h
	${CC} -DMAIN -c -I ${COMMONDIR} -o canonizer_main.${OBJ_EXT} ${COMMONDIR}/canonizer.c
#
smi2mol_main.${OBJ_EXT} : ${COMMONDIR}/smi2mol.c  ${COMMONDIR}/smi2mol.h
	${CC} -DMAIN -c -I ${COMMONDIR} -o smi2mol_main.${OBJ_EXT} ${COMMONDIR}/smi2mol.c
#
# Regulare object modules
#
forio.${OBJ_EXT} : ${COMMONDIR}/forio.c
	${CC} -c -I ${COMMONDIR} -o forio.${OBJ_EXT} ${COMMONDIR}/forio.c
#
canonizer.${OBJ_EXT} : ${COMMONDIR}/canonizer.c  ${COMMONDIR}/canonizer.h
	${CC} -c -I ${COMMONDIR} -o canonizer.${OBJ_EXT} ${COMMONDIR}/canonizer.c
#
matchtest.${OBJ_EXT} : ${COMMONDIR}/matchtest.c
	${CC} -c -I ${COMMONDIR} -o matchtest.${OBJ_EXT} ${COMMONDIR}/matchtest.c
#
local.${OBJ_EXT} : ${COMMONDIR}/local.c ${COMMONDIR}/local.h
	${CC} -c -I ${COMMONDIR} -o local.${OBJ_EXT} ${COMMONDIR}/local.c
#
utilities.${OBJ_EXT} : ${COMMONDIR}/utilities.c  ${COMMONDIR}/utilities.h
	${CC} -c -I ${COMMONDIR} -o utilities.${OBJ_EXT} ${COMMONDIR}/utilities.c
#
smi2mol.${OBJ_EXT} : ${COMMONDIR}/smi2mol.c  ${COMMONDIR}/smi2mol.h
	${CC} -c -I ${COMMONDIR} -o smi2mol.${OBJ_EXT} ${COMMONDIR}/smi2mol.c
#
shortcut.${OBJ_EXT} : ${COMMONDIR}/shortcut.c  ${COMMONDIR}/shortcut.h
	${CC} -c -I ${COMMONDIR} -o shortcut.${OBJ_EXT} ${COMMONDIR}/shortcut.c
#
set.${OBJ_EXT} : ${COMMONDIR}/set.c  ${COMMONDIR}/set.h
	${CC} -c -I ${COMMONDIR} -o set.${OBJ_EXT} ${COMMONDIR}/set.c
#
graph.${OBJ_EXT} : ${COMMONDIR}/graph.c  ${COMMONDIR}/graph.h
	${CC} -c -I ${COMMONDIR} -o graph.${OBJ_EXT} ${COMMONDIR}/graph.c
#
layout.${OBJ_EXT} : ${COMMONDIR}/layout.c  ${COMMONDIR}/layout.h
	${CC} -c -I ${COMMONDIR} -o layout.${OBJ_EXT} ${COMMONDIR}/layout.c
#
denormal.${OBJ_EXT} : ${COMMONDIR}/denormal.c  ${COMMONDIR}/denormal.h
	${CC} -c -I ${COMMONDIR} -o denormal.${OBJ_EXT} ${COMMONDIR}/denormal.c
#
geometry.${OBJ_EXT} : ${COMMONDIR}/geometry.c  ${COMMONDIR}/geometry.h
	${CC} -c -I ${COMMONDIR} -o geometry.${OBJ_EXT} ${COMMONDIR}/geometry.c
#
casutils.${OBJ_EXT} : ${COMMONDIR}/casutils.c  ${COMMONDIR}/casutils.h
	${CC} -c -I ${COMMONDIR} -o casutils.${OBJ_EXT} ${COMMONDIR}/casutils.c
#
pattern.${OBJ_EXT} : ${COMMONDIR}/pattern.c  ${COMMONDIR}/pattern.h
	${CC} -c -I ${COMMONDIR} -o pattern.${OBJ_EXT} ${COMMONDIR}/pattern.c
#
symbol_lists.${OBJ_EXT} : ${COMMONDIR}/symbol_lists.c  ${COMMONDIR}/symbol_lists.h
	${CC} -c -I ${COMMONDIR} -o symbol_lists.${OBJ_EXT} ${COMMONDIR}/symbol_lists.c
#
stereo.${OBJ_EXT} : ${COMMONDIR}/stereo.c  ${COMMONDIR}/stereo.h
	${CC} -c -I ${COMMONDIR} -o stereo.${OBJ_EXT} ${COMMONDIR}/stereo.c
#
perceive.${OBJ_EXT} : ${COMMONDIR}/perceive.c  ${COMMONDIR}/perceive.h
	${CC} -c -I ${COMMONDIR} -o perceive.${OBJ_EXT} ${COMMONDIR}/perceive.c
#
rtutils.${OBJ_EXT} : ${COMMONDIR}/rtutils.c  ${COMMONDIR}/rtutils.h
	${CC} -c -I ${COMMONDIR} -o rtutils.${OBJ_EXT} ${COMMONDIR}/rtutils.c
#
symboltable.${OBJ_EXT} : ${COMMONDIR}/symboltable.c  ${COMMONDIR}/symboltable.h
	${CC} -c -I ${COMMONDIR} -o symboltable.${OBJ_EXT} ${COMMONDIR}/symboltable.c
#
reaccsio.${OBJ_EXT} : ${COMMONDIR}/reaccsio.c  ${COMMONDIR}/reaccsio.h
	${CC} -c -I ${COMMONDIR} -o reaccsio.${OBJ_EXT} ${COMMONDIR}/reaccsio.c
#
didepict.${OBJ_EXT} : ${COMMONDIR}/didepict.c  ${COMMONDIR}/didepict.h
	${CC} -c -I ${COMMONDIR} -o didepict.${OBJ_EXT} ${COMMONDIR}/didepict.c
#
sketch.${OBJ_EXT} : ${COMMONDIR}/sketch.c  ${COMMONDIR}/sketch.h
	${CC} -c -I ${COMMONDIR} -o sketch.${OBJ_EXT} ${COMMONDIR}/sketch.c
#
hashcode.${OBJ_EXT} : ${COMMONDIR}/hashcode.c  ${COMMONDIR}/hashcode.h
	${CC} -c -I ${COMMONDIR} \
       	       -o hashcode.${OBJ_EXT} ${COMMONDIR}/hashcode.c
#
ssmatch.${OBJ_EXT} : ${COMMONDIR}/ssmatch.c  ${COMMONDIR}/ssmatch.h
	${CC} -c -I ${COMMONDIR} \
       	       -o ssmatch.${OBJ_EXT} ${COMMONDIR}/ssmatch.c
#
depictutil.${OBJ_EXT} : ${COMMONDIR}/depictutil.c  ${COMMONDIR}/depictutil.h
	${CC} -c -I ${COMMONDIR} \
       	       -o depictutil.${OBJ_EXT} ${COMMONDIR}/depictutil.c
#
patclean.${OBJ_EXT} : ${COMMONDIR}/patclean.c  ${COMMONDIR}/patclean.h
	${CC} -c -I ${COMMONDIR} \
       	       -o patclean.${OBJ_EXT} ${COMMONDIR}/patclean.c
#
depict.${OBJ_EXT} : ${WINDIR}/depict.c  ${WINDIR}/depict.h
	${CC} -c -I ${WINDIR} \
	       -I ${COMMONDIR} \
       	       -o depict.${OBJ_EXT} ${WINDIR}/depict.c
#
toolkit.${OBJ_EXT} : ${WINDIR}/toolkit.c  ${WINDIR}/toolkit.h
	${CC} -c -I ${WINDIR} \
	       --verbose -I ${COMMONDIR} \
       	       -o toolkit.${OBJ_EXT} ${WINDIR}/toolkit.c
# JNI section
avalon_jni_JNISmi2Mol.${LIB_EXT} : \
	       casutils.${OBJ_EXT} denormal.${OBJ_EXT} depictutil.${OBJ_EXT} forio.${OBJ_EXT} \
	       geometry.${OBJ_EXT} graph.${OBJ_EXT} hashcode.${OBJ_EXT} layout.${OBJ_EXT} \
	       local.${OBJ_EXT} patclean.${OBJ_EXT} perceive.${OBJ_EXT} reaccsio.${OBJ_EXT} \
	       rtutils.${OBJ_EXT} set.${OBJ_EXT} smi2mol.${OBJ_EXT} pattern.${OBJ_EXT} ssmatch.${OBJ_EXT} \
	       stereo.${OBJ_EXT} symbol_lists.${OBJ_EXT} symboltable.${OBJ_EXT} utilities.${OBJ_EXT} \
	       shortcut.${OBJ_EXT} \
	       ${JNIDIR}/avalon_jni_JNISmi2Mol.cpp  ${JNIDIR}/avalon_jni_JNISmi2Mol.h
	${LD} -Wall -D_JNI_IMPLEMENTATION_ ${KILL_AT} \
	      -I ${JDK_INCLUDE} \
              -I ${JDK_MACHINE_INCLUDE} \
	       -I ${COMMONDIR} \
	      -shared ${JNIDIR}/avalon_jni_JNISmi2Mol.cpp \
	       casutils.${OBJ_EXT} denormal.${OBJ_EXT} depictutil.${OBJ_EXT} forio.${OBJ_EXT} \
	       geometry.${OBJ_EXT} graph.${OBJ_EXT} hashcode.${OBJ_EXT} layout.${OBJ_EXT} \
	       local.${OBJ_EXT} patclean.${OBJ_EXT} perceive.${OBJ_EXT} reaccsio.${OBJ_EXT} \
	       rtutils.${OBJ_EXT} set.${OBJ_EXT} smi2mol.${OBJ_EXT} pattern.${OBJ_EXT} ssmatch.${OBJ_EXT} \
	       stereo.${OBJ_EXT} symbol_lists.${OBJ_EXT} symboltable.${OBJ_EXT} utilities.${OBJ_EXT} \
	       shortcut.${OBJ_EXT} \
	      -o avalon_jni_JNISmi2Mol.${LIB_EXT} \
		   ${LD_OPTS}
#
${LIB_PREFIX}JNIMatch.${LIB_EXT} : \
	       canonizer.${OBJ_EXT} casutils.${OBJ_EXT} denormal.${OBJ_EXT} forio.${OBJ_EXT} \
	       geometry.${OBJ_EXT} graph.${OBJ_EXT} hashcode.${OBJ_EXT} layout.${OBJ_EXT} \
	       local.${OBJ_EXT} perceive.${OBJ_EXT} reaccsio.${OBJ_EXT} \
	       rtutils.${OBJ_EXT} set.${OBJ_EXT} smi2mol.${OBJ_EXT} pattern.${OBJ_EXT} ssmatch.${OBJ_EXT} \
	       stereo.${OBJ_EXT} symbol_lists.${OBJ_EXT} symboltable.${OBJ_EXT} utilities.${OBJ_EXT} \
	      ${JNIDIR}/jni_JNIMatch.cpp  ${JNIDIR}/jni_JNIMatch.h
	${LD} -Wall -D_JNI_IMPLEMENTATION_ ${KILL_AT} \
	      -I ${JDK_INCLUDE} \
              -I ${JDK_MACHINE_INCLUDE} \
	       -I ${COMMONDIR} \
              -shared ${JNIDIR}/jni_JNIMatch.cpp \
	       canonizer.${OBJ_EXT} casutils.${OBJ_EXT} denormal.${OBJ_EXT} forio.${OBJ_EXT} \
	       geometry.${OBJ_EXT} graph.${OBJ_EXT} hashcode.${OBJ_EXT} layout.${OBJ_EXT} \
	       local.${OBJ_EXT} perceive.${OBJ_EXT} reaccsio.${OBJ_EXT} \
	       rtutils.${OBJ_EXT} set.${OBJ_EXT} smi2mol.${OBJ_EXT} pattern.${OBJ_EXT} ssmatch.${OBJ_EXT} \
	       stereo.${OBJ_EXT} symbol_lists.${OBJ_EXT} symboltable.${OBJ_EXT} utilities.${OBJ_EXT} \
		   -o ${LIB_PREFIX}JNIMatch.${LIB_EXT} \
		   ${LD_OPTS}
#
JNIWinTools.${LIB_EXT} : \
	       depictutil.${OBJ_EXT} \
	      ${JNIDIR}/jni_JNIWinTools.cpp  ${JNIDIR}/jni_JNIWinTools.h
	${CC} -Wall -mwindows -D_JNI_IMPLEMENTATION_ ${KILL_AT} \
	      -I ${JDK_INCLUDE} \
              -I ${JDK_MACHINE_INCLUDE} \
			  -I ${WINDIR} \
	       -I ${COMMONDIR} \
	       casutils.${OBJ_EXT} \
	       depictutil.${OBJ_EXT} \
	       denormal.${OBJ_EXT} \
	       forio.${OBJ_EXT} \
	       geometry.${OBJ_EXT} \
	       graph.${OBJ_EXT} \
	       local.${OBJ_EXT} \
	       layout.${OBJ_EXT} \
	       perceive.${OBJ_EXT} \
	       reaccsio.${OBJ_EXT} \
	       rtutils.${OBJ_EXT} \
	       set.${OBJ_EXT} \
	       smi2mol.${OBJ_EXT} pattern.${OBJ_EXT} symboltable.${OBJ_EXT} \
	       stereo.${OBJ_EXT} \
	       symbol_lists.${OBJ_EXT} \
	       utilities.${OBJ_EXT} \
		   -shared ${JNIDIR}/jni_JNIWinTools.cpp \
		   -o JNIWinTools.${LIB_EXT} \
		   ${LD_OPTS}
#
ClipboardTools.${LIB_EXT} : \
	      ${JNIDIR}/jni_ClipboardTools.cpp  ${JNIDIR}/jni_ClipboardTools.h
	${CC} -Wall -mwindows -D_JNI_IMPLEMENTATION_ ${KILL_AT} \
	      -I ${JDK_INCLUDE} \
              -I ${JDK_MACHINE_INCLUDE} \
	       -I ${COMMONDIR} \
           -shared ${JNIDIR}/jni_ClipboardTools.cpp \
		   -o ClipboardTools.${LIB_EXT} \
		   ${LD_OPTS}
#
${LIB_PREFIX}JNIDepict.${LIB_EXT} : \
	       casutils.${OBJ_EXT} \
	       denormal.${OBJ_EXT} \
		   depictutil.${OBJ_EXT} \
		   forio.${OBJ_EXT} \
	       geometry.${OBJ_EXT} \
		   graph.${OBJ_EXT} \
		   hashcode.${OBJ_EXT} \
		   layout.${OBJ_EXT} \
	       local.${OBJ_EXT} \
		   patclean.${OBJ_EXT} \
		   perceive.${OBJ_EXT} \
		   reaccsio.${OBJ_EXT} \
		   rtutils.${OBJ_EXT} \
	       set.${OBJ_EXT} \
		   smi2mol.${OBJ_EXT} pattern.${OBJ_EXT} \
		   ssmatch.${OBJ_EXT} \
	       stereo.${OBJ_EXT} \
		   symbol_lists.${OBJ_EXT} \
		   symboltable.${OBJ_EXT} \
		   utilities.${OBJ_EXT} \
	      ${JNIDIR}/jni_JNIDepict.cpp  ${JNIDIR}/jni_JNIDepict.h
	${CC} -Wall -D_JNI_IMPLEMENTATION_ ${KILL_AT} \
	      -I ${JDK_INCLUDE} \
              -I ${JDK_MACHINE_INCLUDE} \
	       -I ${COMMONDIR} \
              -shared ${JNIDIR}/jni_JNIDepict.cpp \
	       casutils.${OBJ_EXT} \
	       denormal.${OBJ_EXT} \
		   depictutil.${OBJ_EXT} \
		   forio.${OBJ_EXT} \
	       geometry.${OBJ_EXT} \
		   graph.${OBJ_EXT} \
		   hashcode.${OBJ_EXT} \
		   layout.${OBJ_EXT} \
	       local.${OBJ_EXT} \
		   patclean.${OBJ_EXT} \
		   perceive.${OBJ_EXT} \
		   reaccsio.${OBJ_EXT} \
		   rtutils.${OBJ_EXT} \
	       set.${OBJ_EXT} \
		   smi2mol.${OBJ_EXT} pattern.${OBJ_EXT} \
		   ssmatch.${OBJ_EXT} \
	       stereo.${OBJ_EXT} \
		   symbol_lists.${OBJ_EXT} \
		   symboltable.${OBJ_EXT} \
		   utilities.${OBJ_EXT} \
		   -o ${LIB_PREFIX}JNIDepict.${LIB_EXT} \
		   ${LD_OPTS}
#
librohde_clib.${LIB_EXT} : \
	aacheck.${OBJ_EXT} \
	fixcharges.${OBJ_EXT} \
	hashcode.${OBJ_EXT} \
	patclean.${OBJ_EXT} \
	ssmatch.${OBJ_EXT} \
	canonizer.${OBJ_EXT} smi2mol.${OBJ_EXT} pattern.${OBJ_EXT} denormal.${OBJ_EXT} casutils.${OBJ_EXT} forio.${OBJ_EXT} geometry.${OBJ_EXT} \
	graph.${OBJ_EXT} layout.${OBJ_EXT} perceive.${OBJ_EXT} \
	reaccsio.${OBJ_EXT} rtutils.${OBJ_EXT} set.${OBJ_EXT} \
	stereo.${OBJ_EXT} symbol_lists.${OBJ_EXT} symboltable.${OBJ_EXT} utilities.${OBJ_EXT} \
	local.${OBJ_EXT}
	${CC} -Wall \
	      -o librohde_clib.${LIB_EXT} \
	      -shared ${COMMONDIR}/struchk.c \
		aacheck.${OBJ_EXT} canonizer.${OBJ_EXT} smi2mol.${OBJ_EXT} pattern.${OBJ_EXT} denormal.${OBJ_EXT} casutils.${OBJ_EXT} fixcharges.${OBJ_EXT} forio.${OBJ_EXT} geometry.${OBJ_EXT} \
		graph.${OBJ_EXT} hashcode.${OBJ_EXT} layout.${OBJ_EXT} \
		patclean.${OBJ_EXT} perceive.${OBJ_EXT} reaccsio.${OBJ_EXT} rtutils.${OBJ_EXT} set.${OBJ_EXT} \
		ssmatch.${OBJ_EXT} stereo.${OBJ_EXT} symbol_lists.${OBJ_EXT} symboltable.${OBJ_EXT} utilities.${OBJ_EXT} \
	    local.${OBJ_EXT} \
	    ${LD_OPTS}
#
struchk${EXE} : struchk_main.${OBJ_EXT} \
	aacheck.${OBJ_EXT} \
	fixcharges.${OBJ_EXT} \
	hashcode.${OBJ_EXT} \
	pattern.${OBJ_EXT} \
	patclean.${OBJ_EXT} \
	ssmatch.${OBJ_EXT} \
	    casutils.${OBJ_EXT} denormal.${OBJ_EXT} forio.${OBJ_EXT} geometry.${OBJ_EXT} \
	    graph.${OBJ_EXT} layout.${OBJ_EXT} perceive.${OBJ_EXT} \
	    reaccsio.${OBJ_EXT} rtutils.${OBJ_EXT} set.${OBJ_EXT} \
	    stereo.${OBJ_EXT} symbol_lists.${OBJ_EXT} symboltable.${OBJ_EXT} utilities.${OBJ_EXT} \
	    local.${OBJ_EXT}
	${LD} -o struchk${EXE} struchk_main.${OBJ_EXT} \
		aacheck.${OBJ_EXT} fixcharges.${OBJ_EXT} forio.${OBJ_EXT} geometry.${OBJ_EXT} \
		graph.${OBJ_EXT} hashcode.${OBJ_EXT} layout.${OBJ_EXT} pattern.${OBJ_EXT} \
		patclean.${OBJ_EXT} perceive.${OBJ_EXT} reaccsio.${OBJ_EXT} set.${OBJ_EXT} \
		ssmatch.${OBJ_EXT} stereo.${OBJ_EXT} symbol_lists.${OBJ_EXT} symboltable.${OBJ_EXT} utilities.${OBJ_EXT} \
	    local.${OBJ_EXT} \
	    ${LD_OPTS}
#
struchk_main.${OBJ_EXT} : ${COMMONDIR}/struchk.c  ${COMMONDIR}/struchk.h
	${CC} -DMAIN -c -I ${COMMONDIR} -o struchk_main.${OBJ_EXT} ${COMMONDIR}/struchk.c
#
aacheck.${OBJ_EXT} : ${COMMONDIR}/aacheck.c
	${CC} -c -I ${COMMONDIR} -o aacheck.${OBJ_EXT} ${COMMONDIR}/aacheck.c
#
fixcharges.${OBJ_EXT} : ${COMMONDIR}/fixcharges.c
	${CC} -c -I ${COMMONDIR} -o fixcharges.${OBJ_EXT} ${COMMONDIR}/fixcharges.c
#
