# Script to remove all accesses to the m_noOffsprings member variable of the
# cell types inside ZFS.
#
# The following cases are identified (see also matchers below):
# - m_cells->a[ ... ].m_noOffsprings
# - cells[ ... ].m_noOffsprings
# - inside_cells->a[ ... ].m_noOffsprings
# - cpu_cells->a[ ... ].m_noOffsprings
#
# The cells/m_cells cases are handled by rewrite_inplace function
# The custom collector cases are handled by rewrite_custom_collector

# Rewrites occurrences inside grid to a_noOffsprings(...), and inside blocks to
# grid().a_noOffsprings(...)
def rewrite_inplace(occurrence, match, ftype):
    if ftype['grid'] or ftype['avg']:
        return "a_noOffsprings(" + occurrence['delimited_str'] + ")"
    else:
        return "grid().a_noOffsprings(" + occurrence['delimited_str'] + ")"

def rewrite_gcell_inplace(occurrence, match, ftype):
    return "parentGId(" + occurrence['delimited_str'] + ")"

# Rewrites ocurrences of custom collectors to a_noOffsprings(collector, ...)
def rewrite_custom_collector(occurrence, match, ftype):
    if ftype['grid'] or ftype['avg']:
        if "cpu_cells" in match['begin_delimiter']:
            return "a_noOffsprings(cpu_cells, " + occurrence['delimited_str'] + ")"
        elif "input_cells" in match['begin_delimiter']:
            return "a_noOffsprings(input_cells, " + occurrence['delimited_str'] + ")"
        else:
            exit("unknown pattern")
    else:
        if "cellsInput" in match['begin_delimiter']:
            return "grid().a_noOffsprings(cellsInput, " + occurrence['delimited_str'] + ")"
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

end_pId_re = "\]\s*\.m_noOffsprings"
end_pId_ggp_re = "\]\s*\.m_noOffsprings\)"
matchers = [
    # matches m_cells->a
    {'begin_delimiter':"m_cells->a\s*\[",
     'end_delimiter':end_pId_re,
     'rewrite_function':rewrite_inplace},
    # matches cells[]
    {'begin_delimiter':"cells\s*\[",
     'end_delimiter':end_pId_re,
     'rewrite_function':rewrite_inplace},
    # # matches cpu_cells
    {'begin_delimiter':"cpu_cells->a\s*\[",
     'end_delimiter':end_pId_re,
     'rewrite_function':rewrite_custom_collector},
    # matches input_cells
    {'begin_delimiter':"input_cells->a\s*\[",
     'end_delimiter':end_pId_re,
     'rewrite_function':rewrite_custom_collector},
    # matches cellsInput
    {'begin_delimiter':"cellsInput\.a\s*\[",
     'end_delimiter':end_pId_re,
     'rewrite_function':rewrite_custom_collector},
    # matches m_cells->a
    {'begin_delimiter':"m_gCells->a\s*\[",
     'end_delimiter':end_pId_re,
     'rewrite_function':rewrite_gcell_inplace},
    # matches cells[]
    {'begin_delimiter':"gCells\s*\[",
     'end_delimiter':end_pId_re,
     'rewrite_function':rewrite_gcell_inplace},
    # matches m_pCells
    {'begin_delimiter':"\*\(m_pCells\s*\[",
     'end_delimiter':end_pId_ggp_re,
     'rewrite_function':rewrite_inplace},
]

tokens = [
            # # matches noOffsprings in grid files (shadowing)
            # {'pattern': r"\bnoOffsprings\b",
            #  'replacement': "noOffsprings_",
            #  'ftype': ['grid'],
            # }
         ]

