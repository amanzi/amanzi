/*
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
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


