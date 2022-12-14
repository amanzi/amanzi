import pytest
import clean

def test_block_comment_bounds1():
    d1 = """
//
// one way of blcok comments
//
//  with slashes
//

and other stuff

//
""".split('\n')

    
    i,j = clean.find_block_comment_bounds(d1, 2)
    print(d1[j-1])
    assert(i == 1)
    assert(j == 6)

def test_block_comment_bounds2():
    d2 = """

/*
 another 
  without slashes
*/

and other stuff

//
""".split('\n')

    i,j = clean.find_block_comment_bounds(d2, 2)
    print(d2[j-1])
    assert(i == 2)
    assert(j == 6)


def test_author1():
    d1a = """
//
// one way of blcok comments
//
// Authors: John Doe
//          Someone Else someon@else.com    
//
// And some trailing stuff
//

and other stuff

//
""".split('\n')

    a_lines, other_lines = clean.find_remove_authors(d1a)
    assert(len(a_lines) == 2)
    assert(a_lines[0] == 'John Doe')
    assert(a_lines[1] == 'Someone Else someon@else.com')
    assert(len(other_lines) == 11)
    print('\n'.join(other_lines))
    assert(other_lines[-4] == "and other stuff")
    


def test_author2():
    d2a = """

/*
 another 
  without slashes

Author: John Doe


*/

and other stuff

//
""".split('\n')

    a_lines, other_lines = clean.find_remove_authors(d2a)
    assert(len(a_lines) == 1)
    assert(a_lines[0] == 'John Doe')
    print('\n'.join(other_lines))
    assert(len(other_lines) == 13)
    assert(other_lines[-4] == "and other stuff")


    
def test_author3():
    d1a = """
//
//
// Authors: John Doe
//          Someone Else someon@else.com    
//
//

and other stuff

//
""".split('\n')

    a_lines, other_lines = clean.find_remove_authors(d1a)
    assert(len(a_lines) == 2)
    assert(a_lines[0] == 'John Doe')
    assert(a_lines[1] == 'Someone Else someon@else.com')
    print('\n'.join(other_lines))
    assert(len(other_lines) == 6)
    assert(other_lines[-4] == "and other stuff")


def test_author4():
    d1a = """
//
//
// Authors: 
//          John Doe
//          Someone Else someon@else.com    
//
//

and other stuff

//
""".split('\n')

    a_lines, other_lines = clean.find_remove_authors(d1a)
    assert(len(a_lines) == 2)
    assert(a_lines[0] == 'John Doe')
    assert(a_lines[1] == 'Someone Else someon@else.com')
    print('\n'.join(other_lines))
    assert(len(other_lines) == 6)
    assert(other_lines[-4] == "and other stuff")
    




def test_copyright1():
    d1a = """
//
// one way of blcok comments
//
// Copyright abc
// is copyright us
// and some others
//
// Authors: John Doe
//          Someone Else someon@else.com    
//
// And some trailing stuff
//

and other stuff

//
""".split('\n')

    a_lines, other_lines = clean.find_remove_copyright(d1a)
    print('\n'.join(a_lines))
    assert(len(a_lines) == 2 + 8)
    assert(a_lines[-4].strip() == 'Authors: John Doe')
    assert(a_lines[-3].strip() == 'Someone Else someon@else.com')
    assert(a_lines[1].strip().startswith('Copyright 2010-202x held jointly'))
    assert(len(other_lines) == 11)
    print('\n'.join(other_lines))
    assert(other_lines[-4] == "and other stuff")
    


def test_copyright2():
    d2a = """

/*
 another 
  without slashes

Copyright abc
is copyright us
and some others

Author: John Doe


*/

and other stuff

//
""".split('\n')

    a_lines, other_lines = clean.find_remove_copyright(d2a)
    assert(len(a_lines) == 1 + 8)
    assert(a_lines[-3].strip() == 'Authors: John Doe')
    print('\n'.join(other_lines))
    assert(len(other_lines) == 13)
    assert(other_lines[-4] == "and other stuff")


    
def test_copyright3():
    d1a = """
//
// License: BSD
//
// Authors: John Doe
//          Someone Else someon@else.com    
//
//

and other stuff

//
""".split('\n')

    a_lines, other_lines = clean.find_remove_copyright(d1a)
    assert(len(a_lines) == 2 + 8)
    assert(a_lines[-4].strip() == 'Authors: John Doe')
    assert(a_lines[-3].strip() == 'Someone Else someon@else.com')
    print('\n'.join(other_lines))
    assert(len(other_lines) == 6)
    assert(other_lines[-4] == "and other stuff")


def test_cornercase_terse():
    d = """
/*
  This is the flow component of the Amanzi code.
  License: BSD
  Authors: Konstantin Lipnikov (version 2) (lipnikov@lanl.gov)
*/

code
""".split('\n')

    a_lines, other_lines = clean.find_remove_copyright(d)
    assert(len(a_lines) == 1 + 8)

    
