# Java class files, use jni directory as a proxy to trigger compilation of all Java sources
# JAVA_LIBS = lib/log4j-1.2.8.jar;lib/commons-logging-1.1.jar;lib/servlet.jar
JAVA_LIBS = ${SUPPORT_LIBS}
#
javaclasses :
	javac -d ../../Binaries/classes -sourcepath . jni/JNIDepict.java
	javac -Xlint:unchecked -d ../../Binaries/classes -sourcepath . jni/JNIMatch.java
	javac -d ../../Binaries/classes -sourcepath . jni/JNIWinTools.java
	javac -d ../../Binaries/classes -sourcepath . jni/UnixTools.java
	javac -d ../../Binaries/classes -sourcepath . avalon/jni/JNISmi2Mol.java
	sed "s/#VERSION_NUMBER#/${VERSION}/" <avalon/jsp/servlets/Mol2ImageServlet.java >avalon/jsp/servlets/Mol2ImageServlet.java_versioned
	cp avalon/jsp/servlets/Mol2ImageServlet.java_versioned avalon/jsp/servlets/Mol2ImageServlet.java
	javac -Xlint:deprecation,unchecked -cp "${JAVA_LIBS}" -d ../../Binaries/classes -sourcepath . avalon/jsp/servlets/Mol2ImageServlet.java
# Create header files for C code
javastubs :
	javah -classpath "../../Binaries/classes" -d ../jni jni.JNIDepict
	javah -classpath "../../Binaries/classes" -d ../jni jni.JNIMatch
	javah -classpath "../../Binaries/classes" -d ../jni jni.JNIWinTools
	javah -classpath "../../Binaries/classes" -d ../jni avalon.jni.JNISmi2Mol
#
