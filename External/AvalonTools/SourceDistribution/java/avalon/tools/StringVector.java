//
//  Copyright (c) 2010, Novartis Institutes for BioMedical Research Inc.
//  All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met: 
//
//     * Redistributions of source code must retain the above copyright 
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following 
//       disclaimer in the documentation and/or other materials provided 
//       with the distribution.
//     * Neither the name of Novartis Institutes for BioMedical Research Inc. 
//       nor the names of its contributors may be used to endorse or promote 
//       products derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
package avalon.tools;

import java.util.Vector;

/**
 * Implements a specialized Vector-like class that can be used to represent
 * sequences of strings, say, in indented XML generation.
 *
 * It also contains some string related utilities related to simple String
 * parsing.
 */
public class StringVector
{
    Vector lines = null;

    public StringVector()
    {
        lines = new Vector();
    }

    /**
     * Utility function to parse a string into an array of its separator
     * delimited components.
     *
     * Treats two consecutive separators as two tokens 
     * (this is different from how StringTokenizer works)
     */
    public static String[] stringToArray(String line, char separator)
    {
        return (new StringVector(line, separator)).toStringArray();
    }

    /**
     * Create a <code>StringVector</code> by parsing a <code>String</code>
     */
    public StringVector(String line, char separator)
    {
        lines = new Vector();
        while (!line.equals(""))
        {
	    int index = line.indexOf(separator);
	    if (index < 0)
	    {
	        lines.add(line);
	        line = "";
	    }
	    else
	    {
	        String token = line.substring(0, index);
	        line = line.substring(line.indexOf(separator)+1);
	        lines.add(token);
	    }
        }
    }

    /**
     * Convert the contents of this <code>StringVector</code> to an
     * array of <code>String</code>s.
     */
    public String[] toStringArray()
    {
        String result[] = new String[lines.size()];
        for (int i=0; i<result.length; i++)
	    result[i] = (String)lines.elementAt(i);
 
        return result;
    }
}
