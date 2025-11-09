# Script to substitute all old MPI function calls by new ones
#
# The following cases are identified (see also matchers below):
# The following will be done:
# - Identify MPI function calls from list ( -> <this file> )
# - If not wrapped, extend the call according to the list ( -> <the other file> )
# - additionally, append a /*wrapped*/ comment (needed that there is no 
#   infinte loop, and it is a good way to indicate that it is already wrapped
# 
#   Before compiling, make sure, that you call
#     sed -ie 's,/\*wrapped\*/,,g'
#   to ensure, that commentented MPI calls do not cause any trouble

#def cropName( name ):
#    if ( name.startswith('\n') ):
#        name = name.split('\n')[1]
#    while ( name.startswith(' ') or name.startswith('&') ):
#        name = name[1:]
#    return name
def cropName( name ):
    while (name.startswith(' ') or name.startswith('&') or name.startswith('\n')):
        name = name[1:]
    return name
 
# Wraps a function, that only need AT_ macro as appendix
def rewrite_name(occurrence, match, ftype):
    begin = match['begin_delimiter'][:-5]+'('
    return begin + occurrence['delimited_str'] + ', AT_ )/*wrapped*/;'

# Wraps a function, that appends AT_ macro aswell as a variable at 
# position pos
def rewrite_nameAndVarname(occurrence, match, ftype):
    begin = match['begin_delimiter'][:-5]+'('
    pos = match['pos']
    content = occurrence['delimited_str'].split(',')
    varname = content[pos] # extract pos-th position of delimited string
    # Cut off potential whitespace or ampersands
    varname = cropName(varname)
    varname = '\"' + str(varname) + '\"'
    return begin + occurrence['delimited_str'] + ', AT_, '+varname+' )/*wrapped*/;'

# Wraps a function, that appends AT_ macro aswell as a some variables at 
# positions pos
def rewrite_nameAndSndRcvName(occurrence, match, ftype):
    begin = match['begin_delimiter'][:-5]+'('
    pos = match['pos']
    content = occurrence['delimited_str'].split(',')
    
    sndname = content[pos[0]] # extract pos-th position of delimited string
    # Cut off potential whitespace or ampersands
    sndname = cropName(sndname) 
    sndname = '\"' + str(sndname) + '\"'
 
    rcvname = content[pos[1]] # extract pos-th position of delimited string
    # Cut off potential whitespace or ampersands
    rcvname = cropName(rcvname) 
    rcvname = '\"' + str(rcvname) + '\"'

    outname = sndname + ', ' + rcvname 
 
    return begin + occurrence['delimited_str'] + ', AT_, '+outname+' )/*wrapped*/;'

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

# Matchers for all routines implemented in zfsmpioverride.cpp, except 
# communicator related stuff, as this is already finished.


endBracket = '\)\s*;'

matchers = [

      # Point-to-point communication
      {'begin_delimiter':'MPI_Send\s*\(',
       'end_delimiter':endBracket,
       'rewrite_function':rewrite_nameAndVarname,
       'pos': 0},
      {'begin_delimiter':'MPI_Isend\s*\(',
       'end_delimiter':endBracket,
       'rewrite_function':rewrite_nameAndVarname,
       'pos': 0},
      {'begin_delimiter':'MPI_Issend\s*\(',
       'end_delimiter':endBracket,
       'rewrite_function':rewrite_nameAndVarname,
       'pos': 0},
      {'begin_delimiter':'MPI_Recv\s*\(',
       'end_delimiter':endBracket,
       'rewrite_function':rewrite_nameAndVarname,
       'pos': 0},
      {'begin_delimiter':'MPI_Irecv\s*\(',
       'end_delimiter':endBracket,
       'rewrite_function':rewrite_nameAndVarname,
       'pos': 0},
      {'begin_delimiter':'MPI_Send_init\s*\(',
       'end_delimiter':endBracket,
       'rewrite_function':rewrite_nameAndVarname,
       'pos': 0},
      {'begin_delimiter':'MPI_Recv_init\s*\(',
       'end_delimiter':endBracket,
       'rewrite_function':rewrite_nameAndVarname,
       'pos': 0},


      # Wait - Test
      {'begin_delimiter':'MPI_Wait\s*\(',
       'end_delimiter':endBracket,
       'rewrite_function':rewrite_name},
      {'begin_delimiter':'MPI_Waitall\s*\(',
       'end_delimiter':endBracket,
       'rewrite_function':rewrite_name},
      {'begin_delimiter':'MPI_Waitsome\s*\(',
       'end_delimiter':endBracket,
       'rewrite_function':rewrite_name},
      {'begin_delimiter':'MPI_Test\s*\(',
       'end_delimiter':endBracket,
       'rewrite_function':rewrite_name},



      # Collectives
      {'begin_delimiter':'MPI_Barrier\s*\(',
       'end_delimiter':'\)\s*;',
       'rewrite_function':rewrite_name},
      {'begin_delimiter':'MPI_Reduce\s*\(',
       'end_delimiter':'\)\s*;',
       'rewrite_function':rewrite_nameAndSndRcvName,
       'pos': [0, 1]},
      {'begin_delimiter':'MPI_Allreduce\s*\(',
       'end_delimiter':'\)\s*;',
       'rewrite_function':rewrite_nameAndSndRcvName,
       'pos': [0, 1]},
      {'begin_delimiter':'MPI_Scatter\s*\(',
       'end_delimiter':'\)\s*;',
       'rewrite_function':rewrite_nameAndSndRcvName,
       'pos': [0, 3]},
      {'begin_delimiter':'MPI_Bcast\s*\(',
       'end_delimiter':endBracket,
       'rewrite_function':rewrite_nameAndVarname,
       'pos': 0},
      {'begin_delimiter':'MPI_Gather\s*\(',
       'end_delimiter':'\)\s*;',
       'rewrite_function':rewrite_nameAndSndRcvName,
       'pos': [0, 3]},
      {'begin_delimiter':'MPI_Gatherv\s*\(',
       'end_delimiter':'\)\s*;',
       'rewrite_function':rewrite_nameAndSndRcvName,
       'pos': [0, 3]},
      {'begin_delimiter':'MPI_Allgather\s*\(',
       'end_delimiter':'\)\s*;',
       'rewrite_function':rewrite_nameAndSndRcvName,
       'pos': [0, 3]},
      {'begin_delimiter':'MPI_Allgatherv\s*\(',
       'end_delimiter':'\)\s*;',
       'rewrite_function':rewrite_nameAndSndRcvName,
       'pos': [0, 3]},
      {'begin_delimiter':'MPI_Alltoall\s*\(',
       'end_delimiter':'\)\s*;',
       'rewrite_function':rewrite_nameAndSndRcvName,
       'pos': [0, 3]},
      {'begin_delimiter':'MPI_Alltoallv\s*\(',
       'end_delimiter':'\)\s*;',
       'rewrite_function':rewrite_nameAndSndRcvName,
       'pos': [0, 4]},
      {'begin_delimiter':'MPI_Exscan\s*\(',
       'end_delimiter':'\)\s*;',
       'rewrite_function':rewrite_nameAndSndRcvName,
       'pos': [0, 1]},



      # MPI Datatypes
      {'begin_delimiter':'MPI_Type_commit\s*\(',
       'end_delimiter':endBracket,
       'rewrite_function':rewrite_name},
      {'begin_delimiter':'MPI_Type_free\s*\(',
       'end_delimiter':endBracket,
       'rewrite_function':rewrite_name},
      {'begin_delimiter':'MPI_Type_create_hindexed\s*\(',
       'end_delimiter':endBracket,
       'rewrite_function':rewrite_name},
      {'begin_delimiter':'MPI_Type_contiguous\s*\(',
       'end_delimiter':endBracket,
       'rewrite_function':rewrite_name},


      # MPI Group
      {'begin_delimiter':'MPI_Group_incl\s*\(',
       'end_delimiter':endBracket,
       'rewrite_function':rewrite_name},
      {'begin_delimiter':'MPI_Group_free\s*\(',
       'end_delimiter':endBracket,
       'rewrite_function':rewrite_name},


      # MISC
      {'begin_delimiter':'MPI_Start\s*\(',
       'end_delimiter':endBracket,
       'rewrite_function':rewrite_name},
      {'begin_delimiter':'MPI_Startall\s*\(',
       'end_delimiter':endBracket,
       'rewrite_function':rewrite_name},
      {'begin_delimiter':'MPI_Get_count\s*\(',
       'end_delimiter':endBracket,
       'rewrite_function':rewrite_name},
      {'begin_delimiter':'MPI_Get_address\s*\(',
       'end_delimiter':endBracket,
       'rewrite_function':rewrite_name},


      # MPI Info
      {'begin_delimiter':'MPI_Info_create\s*\(',
       'end_delimiter':endBracket,
       'rewrite_function':rewrite_name},
      {'begin_delimiter':'MPI_Info_free\s*\(',
       'end_delimiter':endBracket,
       'rewrite_function':rewrite_name},
      {'begin_delimiter':'MPI_Info_get\s*\(',
       'end_delimiter':endBracket,
       'rewrite_function':rewrite_name},
      {'begin_delimiter':'MPI_Info_get_nthkey\s*\(',
       'end_delimiter':endBracket,
       'rewrite_function':rewrite_name},
      {'begin_delimiter':'MPI_Info_get_nkeys\s*\(',
       'end_delimiter':endBracket,
       'rewrite_function':rewrite_name},
      {'begin_delimiter':'MPI_Info_get_valuelen\s*\(',
       'end_delimiter':endBracket,
       'rewrite_function':rewrite_name},


      # MPI File
      {'begin_delimiter':'MPI_File_open\s*\(',
       'end_delimiter':endBracket,
       'rewrite_function':rewrite_name},
      {'begin_delimiter':'MPI_File_delete\s*\(',
       'end_delimiter':endBracket,
       'rewrite_function':rewrite_name},
      {'begin_delimiter':'MPI_File_seek\s*\(',
       'end_delimiter':endBracket,
       'rewrite_function':rewrite_name},
      {'begin_delimiter':'MPI_File_close\s*\(',
       'end_delimiter':endBracket,
       'rewrite_function':rewrite_name},
      {'begin_delimiter':'MPI_File_write_shared\s*\(',
       'end_delimiter':endBracket,
       'rewrite_function':rewrite_name},
      {'begin_delimiter':'MPI_File_iwrite_shared\s*\(',
       'end_delimiter':endBracket,
       'rewrite_function':rewrite_name},


]

tokens = [
            # # matches globalId in grid files (shadowing)
            # {'pattern': r'\bglobalId\b',
            #  'replacement': 'globalId_',
            #  'ftype': ['grid'],
            # }
         ]
