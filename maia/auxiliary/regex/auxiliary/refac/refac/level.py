# Script to remove all accesses to the m_level member variable of the
# cell types inside ZFS.
#
# The following cases are identified (see also matchers below):
# - m_cells->a[ ... ].m_level
# - cells[ ... ].m_level
#
# The cells/m_cells cases are handled by rewrite_inplace function
# The custom collector cases are handled by rewrite_custom_collector

# Rewrites occurrences inside grid to a_level(...), and inside blocks to
# grid().a_level(...)
def rewrite_inplace(occurrence, match, ftype):
    if ftype['grid'] or ftype['avg'] or ftype['cell']:
        return "a_level(" + occurrence['delimited_str'] + ")"
    else:
        return "grid().a_level(" + occurrence['delimited_str'] + ")"

def rewrite_gcell_inplace(occurrence, match, ftype):
    return "a_levelG(" + occurrence['delimited_str'] + ")"

def rewrite_pp_block_inplace(occurrence, match, ftype):
    return "ppblock()->a_level(" + occurrence['delimited_str'] + ")"

def rewrite_pp_grid_inplace(occurrence, match, ftype):
    return "grid->a_level(" + occurrence['delimited_str'] + ")"

# Rewrites ocurrences of custom collectors to a_level(collector, ...)
def rewrite_custom_collector(occurrence, match, ftype):
    if ftype['grid'] or ftype['avg']:
        if "cpu_cells" in match['begin_delimiter']:
            return "a_level(cpu_cells, " + occurrence['delimited_str'] + ")"
        elif "input_cells" in match['begin_delimiter']:
            return "a_level(input_cells, " + occurrence['delimited_str'] + ")"
        else:
            exit("unknown pattern")
    else:
        if "cellsInput" in match['begin_delimiter']:
            return "grid().a_level(cellsInput, " + occurrence['delimited_str'] + ")"
        else:
            exit("unknown pattern")

# Matchers contains a list of things to be found, and the information
# necessary to rewrite them:
# - Either:
#   - begin/end delimiter must be present to select a range:
#     - if begin ends with [ or (, end must begin with ] or ), respectively
#     - otherwise the first occurrence after begin of end will be chosen
#   - token must be present, a single token will be rewritten:
# - a 'rewrite_function':function(matched, file_type) is required
#   - depending on begin/end vs token, a different matched dictionary
#   will be provided
#
# Matched:
# - in the case of begin/end: it provides: begin, delimited, end
#   - if end ends with [ or ( it also provide delimited_end
# - in the case of token: it provides the token

end_level_re = "\]\s*\.m_level\[[^\]]*\]"
matchers = [
    # matches ppblock()->m_cells->a
    {'begin_delimiter':"ppblock()->m_cells->a\s*\[",
     'end_delimiter':end_level_re,
     'rewrite_function':rewrite_pp_block_inplace},
    # matches grid->m_cells->a
    {'begin_delimiter':"grid->m_cells->a\s*\[",
     'end_delimiter':end_level_re,
     'rewrite_function':rewrite_pp_grid_inplace},
    # matches m_cells->a
    {'begin_delimiter':"m_cells->a\s*\[",
     'end_delimiter':end_level_re,
     'rewrite_function':rewrite_inplace},
    # matches cells[]
    {'begin_delimiter':"cells\s*\[",
     'end_delimiter':end_level_re,
     'rewrite_function':rewrite_inplace},
    # # matches cpu_cells
    {'begin_delimiter':"cpu_cells->a\s*\[",
     'end_delimiter':end_level_re,
     'rewrite_function':rewrite_custom_collector},
    # matches input_cells
    {'begin_delimiter':"input_cells->a\s*\[",
     'end_delimiter':end_level_re,
     'rewrite_function':rewrite_custom_collector},
    # matches cellsInput
    {'begin_delimiter':"cellsInput\.a\s*\[",
     'end_delimiter':end_level_re,
     'rewrite_function':rewrite_custom_collector},
    # matches m_gCells->a[]
    {'begin_delimiter':"m_gCells->a\s*\[",
     'end_delimiter':end_level_re,
     'rewrite_function':rewrite_gcell_inplace},
    # matches gCells[]
    {'begin_delimiter':"gCells\s*\[",
     'end_delimiter':end_level_re,
     'rewrite_function':rewrite_gcell_inplace},
]

tokens = [
            # # matches parentId in grid files (shadowing)
            # {'pattern': r"\bparentId\b",
            #  'replacement': "parentId_",
            #  'ftype': ['grid'],
            # }
         ]
