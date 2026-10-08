from format_short_blocks import reformatCppCode

def format_code(snippet: str) -> str:
    lines = snippet.split('\n')
    return '\n'.join(reformatCppCode([line.rstrip() for line in lines]))

def test_if_with_brace_single_line():
    input_code = 'if (true) { do_something(); }'
    expected = 'if (true) {\n  do_something();\n}'
    assert format_code(input_code) == expected

def test_else_with_brace_single_line():
    input_code = 'else { do_something(); }'
    expected = 'else {\n  do_something();\n}'
    assert format_code(input_code) == expected

def test_if_without_brace_short():
    input_code = 'if (true) do_something();'
    expected = 'if (true) do_something();'
    assert format_code(input_code) == expected

def test_if_without_brace_long():
    input_code = (
        'if (some_extremely_long_condition_that_makes_the_line_go_over_100_characters) '
        'do_something_with_long_args();'
    )
    expected = (
        'if (some_extremely_long_condition_that_makes_the_line_go_over_100_characters)\n'
        '  do_something_with_long_args();'
    )
    assert format_code(input_code) == expected

def test_else_without_brace_short():
    input_code = 'else do_something();'
    expected = 'else do_something();'
    assert format_code(input_code) == expected

def test_else_without_brace_long():
    input_code = (
        'else do_something_with_a_really_long_name_that_exceeds_the_line_limit_by_far_with_a_really_long_name_that_exceeds_the_line_limit_by_far();'
    )
    expected = (
        'else\n'
        '  do_something_with_a_really_long_name_that_exceeds_the_line_limit_by_far_with_a_really_long_name_that_exceeds_the_line_limit_by_far();'
    )
    assert format_code(input_code) == expected

def test_for_with_brace():
    input_code = 'for (int i = 0; i < 10; ++i) { do_something(i); }'
    expected = 'for (int i = 0; i < 10; ++i) {\n  do_something(i);\n}'
    assert format_code(input_code) == expected

def test_for_without_brace_short():
    input_code = 'for (int i = 0; i < 3; ++i) do_something(i);'
    expected = 'for (int i = 0; i < 3; ++i) do_something(i);'
    assert format_code(input_code) == expected

def test_if_with_brace_indented():
    input_code = '      if (condition) { do_something(); }'
    expected = '      if (condition) {\n        do_something();\n      }'
    assert format_code(input_code) == expected
    assert False
    
