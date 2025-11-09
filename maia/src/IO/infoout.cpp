// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "infoout.h"

#include <algorithm> // Needed for string replace function
#include <iomanip>   // Needed for formatted output
#include <iostream>  // Needed for ostream
#include "compiler_config.h"
#if defined(MAIA_MS_COMPILER)
#include <Windows.h>
#include <Winsock2.h>
#include <direct.h>
#else
#include <pwd.h>
#include <unistd.h>
#endif
#include <sys/types.h>
#include <time.h> // Needed for generating time strings
#include "COMM/mpioverride.h"
#include "UTIL/functions.h"
#include "environment.h"

// Needed for the different InfoOutFile filetypes
using namespace MAIA_INFOOUT_FILETYPES;
using namespace std;

/** \brief Adds an attribute to the prefix of the XML string.
 * \author Andreas Lintermann
 * \date 13.08.2012
 *
 * \param[in] att The attribute to add, consists of a pair of MStrings.
 * \return The location of the attribute in the vector of pairs.
 */
MInt InfoOut::addAttribute(pair<MString, MString> att) {
  m_buffer->m_prefixAttributes.push_back(att);
  m_buffer->createPrefixMessage();
  return m_buffer->m_prefixAttributes.size() - 1;
}

/** \brief Erases an attribute from the prefix of the XML string.
 * \author Andreas Lintermann
 * \date 13.08.2012
 *
 *  \param[in] attId The ID of the attribute to delete.
 */
void InfoOut::eraseAttribute(MInt attId) {
  m_buffer->m_prefixAttributes.erase(m_buffer->m_prefixAttributes.begin() + attId);
  m_buffer->createPrefixMessage();
}

/** \brief Modifies an attribute of the prefix of the XML string.
 * \author Andreas Lintermann
 * \date 13.08.2012
 *
 *  \param[in] attId The ID of the attribute to modify.
 *  \param[in] att The new attribute to replace the old one, given by a pair of MStrings.
 */
void InfoOut::modifyAttribute(MInt attId, pair<MString, MString> att) {
  m_buffer->m_prefixAttributes[attId] = att;
  m_buffer->createPrefixMessage();
}

/**
 * \brief Generic constructor is used when no information is provided during declaration.
 * \author Michael Schlottke
 * \date June 2012
 */
InfoOut_buffer::InfoOut_buffer()
  : m_rootOnly(false),
    m_domainId(0),
    m_noDomains(1),
    m_minFlushSize(0),
    m_prefixMessage(),
    m_suffixMessage(),
    m_tmpBuffer(),
    m_projectName() {
  // Nothing here
}

/**
 * \brief Parses the string input and returns the string with XML entities escaped
 * \author Michael Schlottke
 * \date April 2012
 * \details This method iterates over each character of the given input string str and replaces relevant XML
 *          entities with their escaped counterparts.
 *          This code is adapted from http://www.mdawson.net/misc/xmlescape.php (Matson Dawson, 2009).
 *
 * \param[in] str Input string that has characters which need escaping.
 * \return Modified string with XML entities escaped.
 */
MString InfoOut_buffer::encodeXml(const MString& inputStr) {
  MChar c;                       // Contains the current character
  ostringstream tmpEncodeBuffer; // Used as a temporary string buffer

  // Create a for loop that uses an iterator to traverse the complete string
  for(MString::const_iterator iter = inputStr.begin(); iter < inputStr.end(); iter++) {
    // Get current character
    c = (MChar)*iter;

    // Use a switch/case statement for the five XML entities
    switch(c) {
      case '"':
        tmpEncodeBuffer << "&quot;";
        break; // Replace double quotes
      case '&':
        tmpEncodeBuffer << "&amp;";
        break; // Replace ampersand
      case '\'':
        tmpEncodeBuffer << "&apos;";
        break; // Replace single quote
      case '<':
        tmpEncodeBuffer << "&lt;";
        break; // Replace less-than sign
      case '>':
        tmpEncodeBuffer << "&gt;";
        break; // Replace greater-than sign
      default:
        tmpEncodeBuffer << c; // By default just append current character
    }
  }

  // Return encoded stream as a string
  return tmpEncodeBuffer.str();
}

/**
 * \brief Creates an XML prefix using the domain id that is prepended to each message.
 * \author Michael Schlottke, Andreas Lintermann
 * \date 13.08.2012
 *
 * Makes use of an array attribute filled by the user to gerenate the XML string.
 *
 */
inline void InfoOut_buffer::createPrefixMessage() {
  // Create temporary stream
  ostringstream tmpStream;

  // Fill stream with formatted domain id
  tmpStream << "<m d=\"" << m_domainId << "\" ";

  for(MInt i = 0; i < (MInt)(m_prefixAttributes.size()); i++)
    tmpStream << m_prefixAttributes[i].first << "=\"" << m_prefixAttributes[i].second << "\" ";

  tmpStream << ">";

  // Set prefix message to tmpBuffer string
  m_prefixMessage = tmpStream.str();

  // Reset buffer
  tmpStream.str("");
}

/**
 * \brief Creates an XML suffix that is appended to each message.
 * \author Michael Schlottke
 * \date April 2012
 */
inline void InfoOut_buffer::createSuffixMessage() { m_suffixMessage = "</m>\n"; }


/**
 * \brief Sets interal state of whether only the root domain (rank 0) should write to file.
 * \author Michael Schlottke
 * \date June 2012
 *
 * \params[in] rootOnly If true, only rank 0 of the specified MPI communicator writes to file.
 * \return The previous internal state (may be stored to return to the previous behavior).
 */
inline MBool InfoOut_buffer::setRootOnly(MBool rootOnly) {
  MBool previousValue = m_rootOnly;
  m_rootOnly = rootOnly;
  return previousValue;
}

/**
 * \brief Sets the minimum buffer length that has to be reached before the buffer is flushed.
 * \author Michael Schlottke
 * \date June 2012
 *
 * \params[in] minFlushSize Minimum buffer length.
 * \return The previous value of the minimum flush size.
 */
inline MInt InfoOut_buffer::setMinFlushSize(MInt minFlushSize) {
  MInt previousValue = m_minFlushSize;
  m_minFlushSize = minFlushSize;
  return previousValue;
}

/**
 * \brief Return an XML header that should written at the beginning of each log file.
 * \author Michael Schlottke
 * \date June 2012
 * \return The XML header.
 */
MString InfoOut_buffer::getXmlHeader() {
  const MInt maxNoChars = 1024;

  // Gets the current hostname
  MChar host[maxNoChars];
  gethostname(host, maxNoChars - 1);
  host[maxNoChars - 1] = '\0';

  // Gets the current username
  MString user;
#if defined(MAIA_MS_COMPILER)
  constexpr MInt INFO_BUFFER_SIZE = 32767;
  TCHAR infoBuf[INFO_BUFFER_SIZE];
  DWORD bufCharCount = INFO_BUFFER_SIZE;
  if(!GetUserName(infoBuf, &bufCharCount)) {
    user = "n/a";
  } else {
    user = infoBuf;
  }
#else
  passwd* p;
  p = getpwuid(getuid());
  if(p) {
    user = MString(p->pw_name);
  } else {
    user = "n/a";
  }
#endif

  // Gets the current directory
  MChar dir[maxNoChars];
#if defined(MAIA_MS_COMPILER)
  _getcwd(dir, maxNoChars - 1);
#else
  if(!getcwd(dir, maxNoChars - 1)) {
    TERM(-1);
  }
#endif
  dir[maxNoChars - 1] = '\0';

  // Gets the current executionCommand
  stringstream executionCommand;
  executionCommand.str("");
  // m_argv and m_argc additionaly described in maia.cpp
#ifndef PVPLUGIN
  executionCommand << Environment::m_argv[0];
  for(MInt n = 1; n < Environment::m_argc; n++) {
    executionCommand << " " << Environment::m_argv[n];
  }
#else
  executionCommand << "paraview plugin was started --> no execution command";
#endif


  // Create start timestamp
  MChar tmpDateTime[128];
  tm* timeInfo;
  time_t rawTime;

  // Get the current time and write it to rawTime
  time(&rawTime);

  // Convert to time struct
  timeInfo = localtime(&rawTime);

  // Format time to string and save to buffer
  strftime(tmpDateTime, 128, "%Y-%m-%d %H:%M:%S", timeInfo);

  // Create temporary buffer
  ostringstream tmpBuffer;

  // Write XML header information to buffer
  tmpBuffer << "<?xml version=\"1.0\" standalone=\"yes\" ?>\n";
  tmpBuffer << "<root>\n";
  tmpBuffer << "<meta name=\"noDomains\" content=\"" << m_noDomains << "\" />\n";
  tmpBuffer << "<meta name=\"dateCreation\" content=\"" << tmpDateTime << "\" />\n";
  tmpBuffer << "<meta name=\"fileFormatVersion\" content=\"" << m_fileFormatVersion << "\" />\n";
  tmpBuffer << "<meta name=\"projectName\" content=\"" << m_projectName << "\" />\n";
  tmpBuffer << "<meta name=\"user\" content=\"" << user << "\" />\n";
  tmpBuffer << "<meta name=\"host\" content=\"" << host << " (" << XSTRINGIFY(MAIA_HOST_STRING) << ")"
            << "\" />\n";
  tmpBuffer << "<meta name=\"dir\" content=\"" << dir << "\" />\n";
  tmpBuffer << "<meta name=\"executionCommand\" content=\"" << executionCommand.str() << "\" />\n";
  tmpBuffer << "<meta name=\"revision\" content=\"" << XSTRINGIFY(MAIA_VERSION_STRING) << "\" />\n";
  tmpBuffer << "<meta name=\"build\" content=\"" << XSTRINGIFY(MAIA_COMPILER_STRING) << " "
            << XSTRINGIFY(MAIA_BUILD_TYPE_STRING) << " (" << MString(XSTRINGIFY(MAIA_COMPILER_VERSION_STRING)) << ")"
            << "\" />\n";


  // Return XML header
  return tmpBuffer.str();
}

/**
 * \brief Return an XML footer that should written at the end of each log file.
 * \author Michael Schlottke
 * \date June 2012
 *
 * \return The XML footer.
 */
MString InfoOut_buffer::getXmlFooter() {
  // Create timestamp
  MChar tmpDateTime[128];
  tm* timeInfo;
  time_t rawTime;

  // Get the current time and write it to rawTime
  time(&rawTime);

  // Convert to time struct
  timeInfo = localtime(&rawTime);

  // Format time to string and save to buffer
  strftime(tmpDateTime, 128, "%Y-%m-%d %H:%M:%S", timeInfo);

  // Create temporary buffer
  ostringstream tmpBuffer;

  // Write XML footer to buffer
  tmpBuffer << "<meta name=\"dateClosing\" content=\"" << tmpDateTime << "\" />\n";
  tmpBuffer << "</root>\n";

  // Return XML footer
  return tmpBuffer.str();
}


/**
 * \brief Generic constructor is used when no information is provided during declaration.
 * \author Michael Schlottke
 * \date June 2012
 */
InfoOut_simpleFileBuffer::InfoOut_simpleFileBuffer()
  : m_isOpen(false), m_rootOnlyHardwired(false), m_filename(), m_file(), m_mpiComm() {
  // Nothing here
}

/**
 * \brief This constructor yields a new instance that can immediately be used to write messages to a regular file.
 * \author Michael Schlottke
 * \date June 2012
 * \details Internally, this constructor just passes the parameters to open() (see open() for more details).
 *
 * \param[in] filename          Filename that should be used for the file.
 * \param[in] mpiComm           MPI communicator which is used to determine rank/domain information.
 * \param[in] rootOnlyHardwired If true, only rank 0 creates a file and writes to it. On all other processors, no
 *                              file is opened and at each flushing of the buffer, the buffer content is discarded.
 * \param[in] projectName       Projectname that should be used for the executed program.
 */
InfoOut_simpleFileBuffer::InfoOut_simpleFileBuffer(const MString& filename, const MString& projectName,
                                                   MPI_Comm mpiComm, MBool rootOnlyHardwired)
  : m_isOpen(false), m_rootOnlyHardwired(false), m_filename(), m_file(), m_mpiComm() {
  open(filename, projectName, mpiComm, rootOnlyHardwired);
}

/**
 * \brief Destructor calls close() to close the file.
 * \author Michael Schlottke
 * \date June 2012
 */
InfoOut_simpleFileBuffer::~InfoOut_simpleFileBuffer() { close(); }

/**
 * \brief Initialization of the file I/O environment.
 * \author Michael Schlottke
 * \date June 2012
 * \details After a successful call to this method the file stream is ready to use. Any previous information
 *          written to the buffer is lost when open is called. This function creates a new file as specified in
 *          filename (deleting any existing files with the same name), and creates the XML prefixes and suffixes
 *          used for each message. It also writes the necessary XML header information to the file.
 *
 * \param[in] filename          Filename that should be used for the file.
 * \param[in] mpiComm           MPI communicator which is used to determine rank/domain information.
 * \param[in] rootOnlyHardwired If true, only rank 0 creates a file and writes to it. On all other processors, no
 *                              file is opened and at each flushing of the buffer, the buffer content is discarded.
 * \param[in] projectName       Projectname that should be used for the executed program.
 */
void InfoOut_simpleFileBuffer::open(const MString& filename, const MString& projectName, MPI_Comm mpiComm,
                                    MBool rootOnlyHardwired) {
  // Open file only if it was not yet done
  if(!m_isOpen) {
    // Set MPI communicator group
    m_mpiComm = mpiComm;

    // Get domain id and number of domains
    MPI_Comm_rank(m_mpiComm, &m_domainId);
    MPI_Comm_size(m_mpiComm, &m_noDomains);

    // Set whether only domain 0 should do any writing (including the creation of a file)
    m_rootOnlyHardwired = rootOnlyHardwired;

    // Only open the file if m_rootOnlyHardwired was not set as true. Otherwise the file state remains closed.
    if(!(m_rootOnlyHardwired && m_domainId != 0)) {
      // Set filename
      m_filename = filename;

      // Set projectName
      m_projectName = projectName;

      // Open file
      m_file.open(m_filename.c_str());

      // Clear internal buffer in order to dismiss any previous input
      str("");

      // Create prefix and suffix messages
      createPrefixMessage();
      createSuffixMessage();

      // Write root and meta information to file
      m_file << getXmlHeader() << flush;

      // Set state variable
      m_isOpen = true;
    }
  }
}

/**
 * \brief Closes the file.
 * \author Michael Schlottke
 * \date June 2012
 * \details Any subsequent write statements to the file stream are discarded after this method is called. After
 *          close() is called, an XML footer is written to the file. Then the file is closed.
 */
void InfoOut_simpleFileBuffer::close(MBool forceClose) {
  // forceClose is not needed here (only kept for interface consistency reasons)
  static_cast<void>(forceClose);

  // Only close file if was opened before
  if(m_isOpen) {
    // Force flushing of the internal buffer
    flushBuffer();

    // Write XML footer to file and flush stream
    m_file << getXmlFooter() << flush;

    // Close file stream
    m_file.close();

    // Set state variable
    m_isOpen = false;
  }
}

/**
 * \brief Flushes the buffer by writing the contents to the file.
 * \author Michael Schlottke
 * \date June 2012
 * \details Sync is called automatically when an "endl" is sent to the stream. At first the buffer content is
 *          wrapped in the prefix and suffix messages, then the entire string is written to the file by calling
 *          flushBuffer(). Finally, the internal buffers are reset.
 *
 * \return Zero by default.
 */
MInt InfoOut_simpleFileBuffer::sync() {
  // Only write if the file was already opened
  if(m_isOpen) {
    // Create formatted string, escape any XML entities in the message, and save to temporary buffer
    m_tmpBuffer << m_prefixMessage << encodeXml(str()) << m_suffixMessage;

    // Only write to file if current buffer length exceeds the minimum size for flushing
    if(m_tmpBuffer.str().length() >= (unsigned)m_minFlushSize) {
      // Write the string to the file and flush the stream
      m_file << m_tmpBuffer.str() << flush;

      // Reset temporary buffer
      m_tmpBuffer.str("");
    }
  }

  // Reset internal buffer
  str("");

  // Default return value for sync()
  return 0;
}

/**
 * \brief Flushes the buffer by writing the contents to the file.
 * \author Michael Schlottke
 * \date June 2012
 * \details Sync is called automatically when an "endl" is sent to the stream. At first the buffer content is
 *          wrapped in the prefix and suffix messages, then the entire string is written to the file. Finally,
 *          the internal buffers are reset.
 *
 * \return Zero by default.
 */
inline void InfoOut_simpleFileBuffer::flushBuffer() {
  // Only write if the file was already opened
  if(m_isOpen) {
    // Write the string to the file and flush the stream
    m_file << m_tmpBuffer.str() << flush;

    // Reset temporary buffer
    m_tmpBuffer.str("");
  }
}


/**
 * \brief Generic constructor is used when no information is provided during declaration.
 * \author Michael Schlottke
 * \date June 2012
 */
InfoOut_mpiFileBuffer::InfoOut_mpiFileBuffer()
  : m_maxMessageLength(0),
    m_mpiWriteBufferSize(0),
    m_mpiWriteBuffer(0),
    m_isOpen(false),
    m_filename(),
    m_mpiComm(),
    m_mpiFileHandle(),
    m_mpiRequest(MPI_REQUEST_NULL) {
  // Nothing here
}

/**
 * \brief This constructor yields a new instance that can immediately be used to write messages to an MPI file.
 * \author Michael Schlottke
 * \date April 2012
 * \details Internally, this constructor just passes the parameters to open() (see open() for more details).
 *
 * \param[in] filename Filename that should be used for the MPI file.
 * \param[in] mpiComm  MPI communicator for which the file should be opened.
 * \param[in] m_projectName Projectname given.
 */
InfoOut_mpiFileBuffer::InfoOut_mpiFileBuffer(const MString& filename, const MString& projectName, MPI_Comm mpiComm)
  : m_maxMessageLength(0),
    m_mpiWriteBufferSize(0),
    m_mpiWriteBuffer(0),
    m_isOpen(false),
    m_filename(),
    m_mpiComm(),
    m_mpiFileHandle(),
    m_mpiRequest(MPI_REQUEST_NULL) {
  open(filename, projectName, mpiComm);
}

/**
 * \brief Destructor calls close() to close the MPI file (if opened) and deletes the MPI write buffer
 * \author Michael Schlottke
 * \date April 2012
 * \details The MPI write buffer #m_mpiWriteBuffer is deleted as well, since it could be that the buffer is
 *          allocated using setMinFlushSize() but never opened with open() (and thus not deleted within close()).
 */
InfoOut_mpiFileBuffer::~InfoOut_mpiFileBuffer() {
  // Close the file
  close();

  // Delete write buffer
  delete[] m_mpiWriteBuffer;
}

/**
 * \brief Initialization of the MPI I/O environment.
 * \author Michael Schlottke
 * \date April 2012
 * \details After a successful call to this method the MPI file stream is ready to use. Any previous information written
 * to the buffer is lost when open is called. This function creates a new file as specified in #filename (deleting any
 * existing files with the same name) and creates the XML prefixes and suffixes used for each message. The root process
 * also writes the necessary XML header information to the file.
 *
 * \param[in] filename Name of the file to open.
 * \param[in] mpiComm MPI communicator for which to open the file.
 * \param[in] projectName Projectname given.
 */
void InfoOut_mpiFileBuffer::open(const MString& filename, const MString& projectName, MPI_Comm mpiComm) {
  // Open file only if it was not yet done
  if(!m_isOpen) {
    // Set MPI communicator group
    m_mpiComm = mpiComm;

    // Get domain id and number of domains
    MPI_Comm_rank(m_mpiComm, &m_domainId);
    MPI_Comm_size(m_mpiComm, &m_noDomains);

    // Set filename
    m_filename = filename;

    // Set projectName
    m_projectName = projectName;

    // Delete log file if one exists
    MPI_File_delete((MChar*)m_filename.c_str(), MPI_INFO_NULL);

    // Set MPI config values
    MInt amode = MPI_MODE_CREATE | MPI_MODE_WRONLY;

    // Open MPI file
    MPI_File_open(m_mpiComm, (MChar*)m_filename.c_str(), amode, MPI_INFO_NULL, &m_mpiFileHandle, AT_);

    // Clear internal buffer in order to dismiss any previous input
    str("");

    // Create prefix and suffix messages and calculate the maximum length for a message
    createPrefixMessage();
    createSuffixMessage();
    m_maxMessageLength = m_maxStringLength - m_prefixMessage.length() - m_suffixMessage.length();

    // Create new write buffer with the following size: length before flushing + max MPI string length
    setMpiWriteBuffer(m_minFlushSize + m_maxStringLength);

    // Write root and meta information to log file (only on MPI root, i.e. domain id == 0)
    if(m_domainId == 0) {
      // Copy buffer string to MPI buffer
      MInt mpiStringLength = getXmlHeader().copy(m_mpiWriteBuffer, m_mpiWriteBufferSize);

      // Write to MPI file
      MPI_File_write_shared(m_mpiFileHandle, m_mpiWriteBuffer, mpiStringLength, MPI_CHAR, MPI_STATUS_IGNORE, AT_);
    }

    // Make sure that initial message is written before file is used
    MPI_Barrier(m_mpiComm, AT_);

    // Set state variable
    m_isOpen = true;
  }
}

/**
 * \brief Closes the MPI file.
 * \author Michael Schlottke
 * \date April 2012
 * \details Any subsequent write statements to the MPI file stream are discarded after this method is called. NOTE: If
 * the MPI file is not closed using this function before MPI_Finalize() is called, an MPI error occurs!
 *
 *          After close() is called, the MPI root process first writes an XML footer to the file. Then the file is
 * closed on all ranks.
 */
void InfoOut_mpiFileBuffer::close(MBool forceClose) {
  // Only perform finalization if this instance was initialized before
  if(m_isOpen) {
    // Only perform cleanup if forceClose is not set to true
    if(!forceClose) {
      // Force flushing of the internal buffer
      flushBuffer();

      // Call MPI_Wait to make sure everything was written to disk
      MPI_Wait(&m_mpiRequest, MPI_STATUS_IGNORE, AT_);

      // Make sure that all processes have finished writing, so that the footer can be written before the file is closed
      MPI_Barrier(m_mpiComm, AT_);

      // Write root information to log file (only on MPI root)
      if(m_domainId == 0) {
        // Copy buffer string to MPI buffer
        MInt mpiStringLength = getXmlFooter().copy(m_mpiWriteBuffer, m_mpiWriteBufferSize);

        // Write to MPI file
        MPI_File_write_shared(m_mpiFileHandle, m_mpiWriteBuffer, mpiStringLength, MPI_CHAR, MPI_STATUS_IGNORE, AT_);
      }

      // Close MPI file
      MPI_File_close(&m_mpiFileHandle, AT_);
    }

    // Delete MPI write buffer
    delete[] m_mpiWriteBuffer;

    // Set pointer to null to avoid problems when the destructor calls delete again
    m_mpiWriteBuffer = 0;

    // Reset state variable
    m_isOpen = false;
  }
}

/**
 * \brief Flushes the buffer if flushing conditions are met.
 * \author Michael Schlottke
 * \date June 2012
 * \details Sync is called automatically when an "endl" is sent to the stream. At first the buffer content is
 *          wrapped in the prefix and suffix messages, then the entire string is written to the MPI file using
 *          nonblocking communication by calling flushBuffer(). Finally the internal buffers are reset.
 *
 * \return Zero by default.
 */
MInt InfoOut_mpiFileBuffer::sync() {
  // Only write to MPI if the MPI file was already opened, and only if we're supposed to on this domain
  if(m_isOpen && !(m_rootOnly && m_domainId != 0)) {
    // Escape XML entities in buffer
    str(encodeXml(str()));

    // Create formatted string and write to temporary buffer
    if(str().length() > (unsigned)m_maxMessageLength) {
      // If buffer content + prefix/suffix message is too long for m_mpiWriteBuffer, use substr and ensure a newline at
      // the end. In substr the '-1' is for the newline character.
      m_tmpBuffer << m_prefixMessage << str().substr(0, m_maxMessageLength - 1) << "\n" << m_suffixMessage;
    } else {
      // Otherwise just copy the whole string
      m_tmpBuffer << m_prefixMessage << str() << m_suffixMessage;
    }

    // Only write to MPI file if current buffer length exceeds the minimum size for flushing
    if(m_tmpBuffer.str().length() >= (unsigned)m_minFlushSize) {
      flushBuffer();
    }
  }

  // Reset internal buffer
  str("");

  // Default return value for sync()
  return 0;
}

/**
 * \brief Flushes the buffer by writing the contents to the MPI file.
 * \author Michael Schlottke
 * \date June 2012
 * \details This method flushes the buffer (i.e. copy the buffer to the MPI write buffer, issue the MPI commands
 *          etc.) if the file is open, and if this processor is supposed to write to the file. It ignores any
 *          settings for the minimum flush size, and assumes that all string formatting was already performed.
 */
inline void InfoOut_mpiFileBuffer::flushBuffer() {
  // Only write to MPI if the MPI file was already opened, and only if we're supposed to on this domain
  if(m_isOpen && !(m_rootOnly && m_domainId != 0)) {
    // Before copying the new string into the write buffer, make sure that the previous MPI write has finished
    MPI_Wait(&m_mpiRequest, MPI_STATUS_IGNORE, AT_);

    // Retrieve actual string length while copying the string from the temporary buffer to the write buffer.
    // The copy command is invoked with the current buffer size to ensure that there is no buffer overflow.
    MInt mpiStringLength = m_tmpBuffer.str().copy(m_mpiWriteBuffer, m_mpiWriteBufferSize);

    // Reset temporary buffer
    m_tmpBuffer.str("");

    // Write mpiStringLength characters to file using nonblocking MPI I/O
    MPI_File_iwrite_shared(m_mpiFileHandle, m_mpiWriteBuffer, mpiStringLength, MPI_CHAR, &m_mpiRequest, AT_);
  }
}

/**
 * \brief Sets the minimum buffer length that has to be reached before the buffer is flushed.
 * \author Michael Schlottke
 * \date June 2012
 * \details This calls InfoOut_buffer::setMinFlushSize and then adds some code to change the size of the MPI
 *          write buffer m_mpiWriteBuffer along with the minFlushSize.
 *
 * \params[in] minFlushSize Minimum buffer length.
 * \return The previous value of the minimum flush size.
 */
inline MInt InfoOut_mpiFileBuffer::setMinFlushSize(MInt minFlushSize) {
  // Call base class method for general stuff
  MInt previousValue = InfoOut_buffer::setMinFlushSize(minFlushSize);

  // Flush current buffer
  flushBuffer();

  // Create new write buffer with the following size: length before flushing + max MPI string length
  setMpiWriteBuffer(m_minFlushSize + m_maxStringLength);

  return previousValue;
}

/**
 * \brief Delete the current MPI write buffer and allocate a new one.
 * \author Michael Schlottke
 * \date June 2012
 * \details This method calls MPI_Wait to make sure that all write statements have completed before the buffer is
 *          deleted. Then m_mpiWriteBuffer is deleted and afterwards recreated with the new buffer size.
 *
 * \params[in] newBufferSize The new size of the write buffer
 * \return Previous buffer size
 */
MInt InfoOut_mpiFileBuffer::setMpiWriteBuffer(MInt newBufferSize) {
  // Make sure that the write buffer is empty
  MPI_Wait(&m_mpiRequest, MPI_STATUS_IGNORE, AT_);

  // Delete old write buffer
  delete[] m_mpiWriteBuffer;

  // Create new write buffer
  m_mpiWriteBuffer = new MChar[newBufferSize];

  // Save old buffer size to temporary variable and new buffer size to member variable
  MInt oldBufferSize = m_mpiWriteBufferSize;
  m_mpiWriteBufferSize = newBufferSize;

  // Return the previous buffer size
  return oldBufferSize;
}


/**
 * \brief Default constructor does basically nothing.
 * \author Michael Schlottke
 * \date April 2012
 */
InfoOut_streamBuffer::InfoOut_streamBuffer()
  : m_isInitialized(false), m_printDomainId(false), m_output(0), m_mpiComm(), m_isDisabled(false) {
  // Nothing here
}

/**
 * \brief Constructor passes the parameters to initialize().
 * \author Michael Schlottke
 * \date April 2012
 * \details After an object is instantiated using this constructor, the buffer is ready to use. For detailed information
 * on what is done, please refer to initialize(). \param[in] os Output stream that should be used for output. \param[in]
 * mpiComm MPI communicator for which the domain information should be gathered. \param[in] printDomainId Determines
 * whether the domain id should be prepended to each message.
 */
InfoOut_streamBuffer::InfoOut_streamBuffer(ostream* os, MPI_Comm mpiComm, MBool printDomainId)
  : m_isInitialized(false), m_printDomainId(false), m_output(0), m_mpiComm(), m_isDisabled(false) {
  initialize(os, mpiComm, printDomainId);
}

/**
 * \brief Initializes the buffer to make it ready for use.
 * \author Michael Schlottke
 * \date April 2012
 * \details The output stream that will be written to is assigned and the domain id determined. If specified,
 *          a prefix is generated containing the domain id. This method may be called multiple times to change
 *          the desired behavior of this buffer.
 *
 * \param[in] os Output stream that should be used for output.
 * \param[in] mpiComm MPI communicator for which the domain information should be gathered.
 * \param[in] printDomainId Determines whether the domain id should be prepended to each message.
 */
void InfoOut_streamBuffer::initialize(ostream* os, MPI_Comm mpiComm, MBool printDomainId) {
  // Delete any previous (non-empty) strings stored in the buffer that were not flushed yet by calling sync
  if(!str().empty()) {
    sync();
  }

  // Assign ostream to m_output
  m_output = os;

  // Set MPI communicator group
  m_mpiComm = mpiComm;

  // Save whether the domain id should be prepended to each message
  m_printDomainId = printDomainId;

  // Get domain id
  MPI_Comm_rank(m_mpiComm, &m_domainId);

  // Create prefix string including the domain id if specified
  if(m_printDomainId) {
    // Create temporary stream
    ostringstream tmpStream;

    // Fill tmpStream with formatted domain id
    tmpStream << "[" << setfill('0') << setw(7) << m_domainId << "] ";

    // Set prefix message to tmpStream string
    m_prefixMessage = tmpStream.str();
  } else {
    // Otherwise use an empty prefix message
    m_prefixMessage = "";
  }

  // Set state variable
  m_isInitialized = true;

  // Set isDisabled
  m_isDisabled = false;
}

/**
 * \Writes the domain_Id in front of every  line.
 * \If \n is written in the middle of a stream it will be replaced with the domain_Id
 * \and written in the following line so you have the expected output.
 * \date November 2012
 */

MString InfoOut_streamBuffer::addPrefix(const MString& inputStr) {
  MChar c;
  ostringstream tmpEncodeStreambuffer;

  for(MString::const_iterator iter = inputStr.begin(); iter < inputStr.end(); iter++) {
    c = (MChar)*iter;
    tmpEncodeStreambuffer << c;

    if(c == '\n' && iter != inputStr.end() - 1) {
      tmpEncodeStreambuffer << m_prefixMessage;
    }
  }

  return tmpEncodeStreambuffer.str();
}


/**
 * \brief Flushes the buffer by writing the contents to the asssociated output stream.
 * \author Michael Schlottke
 * \date April 2012
 * \details Sync is called automatically when an "endl" is sent to the stream. It preprends each message with the
 *          optional message prefix and writes the contents of the buffer to the output stream, which is then flushed
 *          itself. Finally, all internal buffers are reset.
 *
 * \return Zero by default.
 */
MInt InfoOut_streamBuffer::sync() {
  // Only write actual output if streamBuffer is initialized, and if we are supposed to on this processor
  if(m_isInitialized && !(m_rootOnly && m_domainId != 0) && !m_isDisabled) {
    // Create formatted string
    m_tmpBuffer << m_prefixMessage << addPrefix(str());

    // Write temporary buffer to output stream and flush it
    *m_output << m_tmpBuffer.str() << flush;
  }

  // Reset stringbuf of InfoOut_streamBuffer
  str("");

  // Reset temporary buffer
  m_tmpBuffer.str("");

  // Default return value for sync()
  return 0;
}

/**
 * \brief This function does nothing for this class, since as of now there is no intermediate buffering
 * \author Michael Schlottke
 * \date June 2012
 */
inline void InfoOut_streamBuffer::flushBuffer() {
  // Nothing here
}


/**
 * \brief Default constructor creates (virtual) file that cannot yet be used.
 * \author Michael Schlottke
 * \date April 2012
 * \details When this constructor is used, open() needs to be called before the stream can be used.
 */
InfoOutFile::InfoOutFile() : m_fileType(0), m_isOpen(false) {
  // Nothing here
}

/**
 * \brief Constructor creates InfoOut_mpiFileBuffer buffer and calls ostream constructor with reference to it.
 * \author Michael Schlottke
 * \date April 2012
 * \details When this constructor is used, the stream is immediately ready to use. For information about the paramters,
 *          please have a look at InfoOut_mpiFileBuffer::open.
 *
 * \param[in] filename Name of the file to open.
 * \param[in] mpiComm MPI communicator for which to open the file.
 * \param[in] projectName Projectname given.
 */
InfoOutFile::InfoOutFile(const MString& filename, const MString& projectName, MInt fileType, MPI_Comm mpiComm,
                         MBool rootOnlyHardwired)
  : m_fileType(0), m_isOpen(false) {
  open(filename, projectName, fileType, mpiComm, rootOnlyHardwired);
}

/**
 * \brief Destructor closes the stream.
 * \author Michael Schlottke
 * \date June 2012
 */
InfoOutFile::~InfoOutFile() { close(); }

/**
 * \brief Opens a file by passing the parameters to InfoOut_<xyz>FileBuffer::open(...).
 * \author Michael Schlottke
 * \date April 2012
 * \details The parameter fileType can be any of the constants defined in namespace MAIA_INFOOUT_FILETYPES. This
 *          method then creates a new internal buffer anbd passes along the parameters.
 *
 * \param[in] filename Name of the file to open.
 * \param[in] fileType Type of file that should be opened. Can be any of the constants defined in
 *                     MAIA_INFOOUT_FILETYPES.
 * \param[in] mpiComm MPI communicator for which to open the file.
 * \param[in] rootOnlyHardwired If true, only rank 0 creates a file and writes to it. On all other processors, no
 *                              file is opened and at each flushing of the buffer, the buffer content is discarded.
 *                              This parameter makes only sense when using MAIA_INFOOUT_SIMPLE_FILE and is ignored
 *                              otherwise.
 * \param[in] projectName Name of the executed program.
 */
void InfoOutFile::open(const MString& filename, const MString& projectName, MInt fileType, MPI_Comm mpiComm,
                       MBool rootOnlyHardwired) {
  // Only open file if it was not yet opened
  if(!m_isOpen) {
    // Save file type to member variable
    m_fileType = fileType;

    // Create a new buffer object depending on the specified file type
    switch(m_fileType) {
      case MAIA_INFOOUT_SIMPLE_FILE:
        // Open a simple file
        m_buffer = new InfoOut_simpleFileBuffer(filename, projectName, mpiComm, rootOnlyHardwired);
        break;
      case MAIA_INFOOUT_MPI_FILE:
        // Open an MPI file
      default:
        // This is also the default behavior if no valid fileType was specified (or none at all).
        m_fileType = MAIA_INFOOUT_MPI_FILE;
        m_buffer = new InfoOut_mpiFileBuffer(filename, projectName, mpiComm);
        break;
    }

    // Associate the stream with the newly created buffer
    rdbuf(m_buffer);

    // Set state variable
    m_isOpen = true;
  }
}

/**
 * \brief Pass the close call to the respective internal buffer.
 * \author Michael Schlottke
 * \date June 2012
 * \details All attempts to write to the stream after closing it will be discarded.
 */
void InfoOutFile::close(MBool forceClose) {
  // Only close file if was already opened
  if(m_isOpen) {
    // Determine correct cast by evaluating internal file type variable
    switch(m_fileType) {
      case MAIA_INFOOUT_SIMPLE_FILE:
        static_cast<InfoOut_simpleFileBuffer*>(m_buffer)->close(forceClose);
        break;
      case MAIA_INFOOUT_MPI_FILE:
        static_cast<InfoOut_mpiFileBuffer*>(m_buffer)->close(forceClose);
        break;
      default: {
        mTerm(1, AT_, "Unknown file type");
      }
    }

    // Delete internal buffer to prevent memory leaks
    delete m_buffer;

    // Set state variable
    m_isOpen = false;
  }
}

/**
 * \brief Sets interal state of whether only the root domain (rank 0) should write to file.
 * \author Michael Schlottke
 * \date June 2012
 *
 * \params[in] rootOnly If true, only rank 0 of the specified MPI communicator writes to file.
 * \return The previous internal state (may be stored to return to the previous behavior).
 */
MBool InfoOutFile::setRootOnly(MBool rootOnly) { return m_buffer->setRootOnly(rootOnly); }

/**
 * \brief Sets the minimum buffer length that has to be reached before the buffer is flushed.
 * \author Michael Schlottke
 * \date June 2012
 * \details Flushing the buffer means that the contents of the buffer in memory is written to the file. If the
 *          file stream was not opened yet, this method just returns 0 and does nothing else.
 *
 * \params[in] minFlushSize Minimum buffer length.
 * \return The previous value of the minimum flush size.
 */
MInt InfoOutFile::setMinFlushSize(MInt minFlushSize) {
  if(m_isOpen) {
    return m_buffer->setMinFlushSize(minFlushSize);
  } else {
    return 0;
  }
}


/**
 * \brief Default construtor creates InfoOutStream that is not (yet) usable.
 * \author Michael Schlottke
 * \date April 2012
 * \details When this constructor is used, initialize() needs to be called before the stream can be used.
 */
InfoOutStream::InfoOutStream() : m_isInitialized(false) {
  // Nothing here
}

/**
 * \brief Constructor creates InfoOut_streamBuffer buffer and calls ostream constructor with reference to it.
 * \author Michael Schlottke
 * \date April 2012
 * \details When this constructor is used, the stream is immediately ready to use. For information about the paramters,
 *          please have a look at InfoOut_streamBuffer::initialize.
 *
 * \param[in] os Output stream that should be used for output.
 * \param[in] mpiComm MPI communicator for which the domain information should be gathered.
 * \param[in] printDomainId Determines whether the domain id should be prepended to each message.
 */
InfoOutStream::InfoOutStream(ostream* os, MPI_Comm mpiComm, MBool printDomainId) : m_isInitialized(false) {
  initialize(os, mpiComm, printDomainId);
}


/**
 * \brief Destructor deletes the internal buffer object.
 * \author Michael Schlottke
 * \date June 2012
 */
InfoOutStream::~InfoOutStream() {
  if(m_isInitialized) {
    delete m_buffer;
    m_isInitialized = false;
  }
}

/**
 * \brief Passes the parameters to a call of InfoOut_streamBuffer::initialize().
 * \author Michael Schlottke
 * \date April 2012
 *
 * \param[in] os Output stream that should be used for output.
 * \param[in] mpiComm MPI communicator for which the domain information should be gathered.
 * \param[in] printDomainId Determines whether the domain id should be prepended to each message.
 */
// Pass call to initialize to streamBuffer
void InfoOutStream::initialize(ostream* os, MPI_Comm mpiComm, MBool printDomainId) {
  if(!m_isInitialized) {
    m_buffer = new InfoOut_streamBuffer(os, mpiComm, printDomainId);

    // Associate the stream with the newly created buffer
    rdbuf(m_buffer);

    // Set state variable
    m_isInitialized = true;
  }
}

/**
 * \brief Sets interal state of whether only the root domain (rank 0) should write to file.
 * Author Michael Schlottke
 * \date June 2012
 *
 * \params[in] rootOnly If true, only rank 0 of the specified MPI communicator writes to file.
 * \return The previous internal state (may be stored to return to the previous behavior).
 */
MBool InfoOutStream::setRootOnly(MBool rootOnly) { return m_buffer->setRootOnly(rootOnly); }
