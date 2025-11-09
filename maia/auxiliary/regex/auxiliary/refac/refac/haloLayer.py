# Script to remove all accesses to the m_haloLayer member variable of the
# cell types inside ZFS.
#
# The following cases are identified (see also matchers below):
# - m_cells->a[ ... ].m_haloLayer
# - cells[ ... ].m_haloLayer
#
# The cells/m_cells cases are handled by rewrite_inplace function
# The custom collector cases are handled by rewrite_custom_collector

# Rewrites occurrences inside grid to a_haloLayer(...), and inside blocks to
# grid().a_haloLayer(...)
def rewrite_inplace(occurrence, match, ftype):
    if ftype['grid'] or ftype['avg']:
        return "a_haloLayer(" + occurrence['delimited_str'] + ")"
    else:
        return "grid().a_haloLayer(" + occurrence['delimited_str'] + ")"

def rewrite_gcell_inplace(occurrence, match, ftype):
    return "a_haloLayerG(" + occurrence['delimited_str'] + ")"

def rewrite_particle_inplace(occurrence, match, ftype):
    if "m_blockPtr" in match['begin_delimiter']:
        return "m_blockPtr->grid().a_haloLayer(" + occurrence['delimited_str'] + ")"

def rewrite_gcell_inplace(occurrence, match, ftype):
    return "a_haloLayerG(" + occurrence['delimited_str'] + ")"

def rewrite_pp_block_inplace(occurrence, match, ftype):
    return "ppblock()->a_haloLayer(" + occurrence['delimited_str'] + ")"

def rewrite_pp_grid_inplace(occurrence, match, ftype):
    return "grid->a_haloLayer(" + occurrence['delimited_str'] + ")"

# Rewrites ocurrences of custom collectors to parentId(collector, ...)
def rewrite_custom_collector(occurrence, match, ftype):
    if ftype['grid'] or ftype['avg']:
        if "cpu_cells" in match['begin_delimiter']:
            return "a_haloLayer(cpu_cells, " + occurrence['delimited_str'] + ")"
        elif "input_cells" in match['begin_delimiter']:
            return "a_haloLayer(input_cells, " + occurrence['delimited_str'] + ")"
        else:
            exit("unknown pattern")
    else:
        if "cellsInput" in match['begin_delimiter']:
            return "grid().a_haloLayer(cellsInput, " + occurrence['delimited_str'] + ")"
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

end_haloLayer_re = "\]\s*\.m_haloLayer\s*\[\s*0\s*\]"
end_haloLayer_re2 = "\]\s*\.m_haloLayer\s*\[\s*j\s*\]"
end_haloLayer_re3 = "\]\s*\.m_haloLayer\s*\[\s*i\s*\]"
end_haloLayer_re4 = "\]\s*\.m_haloLayer\s*\[\s*dim\s*\]"
end_haloLayer_re5 = "\]\s*\.m_haloLayer\s*\[\s*spaceId\s*\]"
end_haloLayer_re6 = "\]\s*\.m_haloLayer\s*\[\s*k\s*\]"
end_haloLayer_re7 = "\]\s*\.m_haloLayer\s*\[\s*1\s*\]"
end_haloLayer_re8 = "\]\s*\.m_haloLayer\s*\[\s*d\s*\]"
end_haloLayer_re9 = "\]\s*\.m_haloLayer\s*\[\s*spaceID\s*\]"
end_haloLayer_re10 = "\]\s*\.m_haloLayer\s*\[\s*spaceId1\s*\]"
end_haloLayer_re11 = "\]\s*\.m_haloLayer\s*\[\s*2\s*\]"
matchers = [
    # matches ppblock()->m_cells->a
    {'begin_delimiter':"ppblock()->m_cells->a\s*\[",
     'end_delimiter':end_haloLayer_re,
     'rewrite_function':rewrite_pp_block_inplace},
    # matches grid->m_cells->a
    {'begin_delimiter':"grid->m_cells->a\s*\[",
     'end_delimiter':end_haloLayer_re,
     'rewrite_function':rewrite_pp_grid_inplace},
    # matches m_cells->a
    {'begin_delimiter':"m_cells->a\s*\[",
     'end_delimiter':end_haloLayer_re,
     'rewrite_function':rewrite_inplace},
    # matches cells[]
    {'begin_delimiter':"cells\s*\[",
     'end_delimiter':end_haloLayer_re,
     'rewrite_function':rewrite_inplace},
    # matches m_cells->a
    {'begin_delimiter':"m_cells->a\s*\[",
     'end_delimiter':end_haloLayer_re2,
     'rewrite_function':rewrite_inplace},
    {'begin_delimiter':"cells\s*\[",
     'end_delimiter':end_haloLayer_re2,
     'rewrite_function':rewrite_inplace},
    # matches m_cells->a
    {'begin_delimiter':"m_cells->a\s*\[",
     'end_delimiter':end_haloLayer_re3,
     'rewrite_function':rewrite_inplace},
    # matches m_cells->a
    {'begin_delimiter':"cells\s*\[",
     'end_delimiter':end_haloLayer_re3,
     'rewrite_function':rewrite_inplace},
    # matches m_cells->a
    {'begin_delimiter':"m_cells->a\s*\[",
     'end_delimiter':end_haloLayer_re4,
     'rewrite_function':rewrite_inplace},
    # matches m_cells->a
    {'begin_delimiter':"cells\s*\[",
     'end_delimiter':end_haloLayer_re4,
     'rewrite_function':rewrite_inplace},
    # matches m_cells->a
    {'begin_delimiter':"cells\s*\[",
     'end_delimiter':end_haloLayer_re5,
     'rewrite_function':rewrite_inplace},
    # matches m_cells->a
    {'begin_delimiter':"cells\s*\[",
     'end_delimiter':end_haloLayer_re6,
     'rewrite_function':rewrite_inplace},
    # matches m_cells->a
    {'begin_delimiter':"m_cells->a\s*\[",
     'end_delimiter':end_haloLayer_re6,
     'rewrite_function':rewrite_inplace},
    # matches m_cells->a
    {'begin_delimiter':"cells\s*\[",
     'end_delimiter':end_haloLayer_re7,
     'rewrite_function':rewrite_inplace},
    # matches m_cells->a
    {'begin_delimiter':"m_cells->a\s*\[",
     'end_delimiter':end_haloLayer_re7,
     'rewrite_function':rewrite_inplace},
    # matches m_cells->a
    {'begin_delimiter':"cellsInput.a\s*\[",
     'end_delimiter':end_haloLayer_re2,
     'rewrite_function':rewrite_custom_collector},
    # matches m_cells->a
    {'begin_delimiter':"\*\(m_blockPtr->cells\s*\[",
     'end_delimiter':"\]\s*\.m_haloLayer\s*\)",
     'rewrite_function':rewrite_particle_inplace},
    # matches m_pCells
    {'begin_delimiter':"\*\(m_pCells\s*\[",
     'end_delimiter':"\]\s*\.m_haloLayer\s*\)",
     'rewrite_function':rewrite_inplace},
    {'begin_delimiter':"cells\s*\[",
     'end_delimiter':end_haloLayer_re8,
     'rewrite_function':rewrite_inplace},
    {'begin_delimiter':"cells\s*\[",
     'end_delimiter':end_haloLayer_re9,
     'rewrite_function':rewrite_inplace},
    {'begin_delimiter':"cells\s*\[",
     'end_delimiter':end_haloLayer_re10,
     'rewrite_function':rewrite_inplace},
    {'begin_delimiter':"m_cells->a\s*\[",
     'end_delimiter':end_haloLayer_re11,
     'rewrite_function':rewrite_inplace},
    {'begin_delimiter':"cpu_cells->a\s*\[",
     'end_delimiter':end_haloLayer_re2,
     'rewrite_function':rewrite_custom_collector},
    # matches input_cells
    {'begin_delimiter':"input_cells->a\s*\[",
     'end_delimiter':end_haloLayer_re2,
     'rewrite_function':rewrite_custom_collector},
    # matches m_cells->a
    {'begin_delimiter':"m_gCells->a\s*\[",
     'end_delimiter':end_haloLayer_re,
     'rewrite_function':rewrite_gcell_inplace},
    # matches cells[]
    {'begin_delimiter':"gCells\s*\[",
     'end_delimiter':end_haloLayer_re,
     'rewrite_function':rewrite_gcell_inplace},
    # matches m_cells->a
    # {'begin_delimiter':"m_gLevel\s*\[",
    #  'end_delimiter':"\]\s*",
    #  'rewrite_function':rewrite_gcell_inplace}
]

tokens = [
            # # matches parentId in grid files (shadowing)
            # {'pattern': r"\bparentId\b",
            #  'replacement': "parentId_",
            #  'ftype': ['grid'],
            # }
         ]
