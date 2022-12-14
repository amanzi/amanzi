/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

/*
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

/*
 * $Id: DOMTreeErrorReporter.cpp 471735 2006-11-06 13:53:58Z amassari $
 */

// ---------------------------------------------------------------------------
//  Includes
// ---------------------------------------------------------------------------
#include <xercesc/dom/DOMParseException.hpp>
#include "DOMTreeErrorReporter.hpp"
#include <iostream>
#include <stdlib.h>
#include <string.h>


void DOMTreeErrorReporter::warning(const DOMParseException&)
{
    //
    // Ignore all warnings.
    //
}

void DOMTreeErrorReporter::error(const DOMParseException& toCatch)
{
#define XSTR(s) STR(s)
#define STR(s) #s

    fSawErrors = true;
    std::cerr << "Error at file \"" << XSTR(toCatch.getSystemId())
		 << "\", line EIB " << toCatch.getLineNumber()
		 << ", column " << toCatch.getColumnNumber()
         << "\n   Message: " << XSTR(toCatch.getMessage()) << std::endl;
}

void DOMTreeErrorReporter::fatalError(const DOMParseException& toCatch)
{
#define XSTR(s) STR(s)
#define STR(s) #s

    fSawErrors = true;
    std::cerr << "Fatal Error at file \"" << XSTR(toCatch.getSystemId())
		 << "\", line EIB " << toCatch.getLineNumber()
		 << ", column " << toCatch.getColumnNumber()
         << "\n   Message: " << XSTR(toCatch.getMessage()) << std::endl;
}

void DOMTreeErrorReporter::resetErrors()
{
    fSawErrors = false;
}
