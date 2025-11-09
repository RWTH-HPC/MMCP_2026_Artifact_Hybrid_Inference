# Script to remove all accesses to the m_delete member variable of the
# cell types inside ZFS.
#
# The following cases are identified (see also matchers below):
# - m_cells->a[ ... ].m_delete
# - cells[ ... ].m_delete
# - inside_cells->a[ ... ].m_delete
# - cpu_cells->a[ ... ].m_delete
#
# The cells/m_cells cases are handled by rewrite_inplace function
# The custom collector cases are handled by rewrite_custom_collector

# Rewrites occurrences inside blocks to delete(...)
def rewrite_inplace(occurrence, match, ftype):
    if ftype['bnd']:
        if "block" in match['begin_delimiter']:
            return "m_block->a_isToDelete(" + occurrence['delimited_str'] + ")"
        else: 
            return "m_block->a_isToDelete(" + occurrence['delimited_str'] + ")"
    else:
        return "a_isToDelete(" + occurrence['delimited_str'] + ") "


# Rewrites ocurrences of custom collectors to delete(collector, ...)
def rewrite_custom_collector(occurrence, match, ftype):
    if ftype['grid'] or ftype['avg']:
        if "cpu_cells" in match['begin_delimiter']:
            return "delete(cpu_cells, " + occurrence['delimited_str'] + ")"
        elif "input_cells" in match['begin_delimiter']:
            return "a_isToDelete(input_cells, " + occurrence['delimited_str'] + ")"
        else:
            exit("unknown pattern")
    else:
        if "cellsInput" in match['begin_delimiter']:
            return "grid().a_isToDelete(cellsInput, " + occurrence['delimited_str'] + ")"
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

end_pId_re = "\]\s*\.m_delete\\b"
end_pId_ggp_re = "\]\s*\.m_delete\)"
matchers = [
     # matches m_block->m_cells->a
    {'begin_delimiter':"m_block->m_cells->a\s*\[",
     'end_delimiter':end_pId_re,
     'rewrite_function':rewrite_inplace},
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
]

tokens = [
            # # matches delete in grid files (shadowing)
            # {'pattern': r"\bdelete\b",
            #  'replacement': "delete_",
            #  'ftype': ['grid'],
            # }
         ]
