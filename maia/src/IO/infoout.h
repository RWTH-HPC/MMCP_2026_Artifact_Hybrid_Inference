// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef INFOOUT_H
#define INFOOUT_H

#include <fstream> // Needed for ofstream
#include <sstream> // Needed for ostringstream
#include <vector>
#include "COMM/mpioverride.h" // Needed for MPI functionality
#include "INCLUDE/maiatypes.h"
#include "COMM/globalmpiinfo.h"

#ifdef _SX
#include <sys/socket.h>
#endif


/**
 * \brief Namespace to hold all supported filetypes within the InfoOutFile
 * \author Michael Schlottke
 * \date June 2012
 */
namespace MAIA_INFOOUT_FILETYPES {
const MInt MAIA_INFOOUT_MPI_FILE = 0;    //!< Use a single file for all domains (MPI I/O)
const MInt MAIA_INFOOUT_SIMPLE_FILE = 1; //!< Use a physical file for each domain
} // namespace MAIA_INFOOUT_FILETYPES

/**
 * \brief Base class for all InfoOut_<xyz>Buffer classes.
 * \author Michael Schlottke
 * \date June 2012
 * \details This class holds all relevant member variables and functions that are needed for all buffer types, i.e.
 *          formattting routines, variables to store domain information etc. Basically everything that is needed in
 *          all (or at least most) specialized buffer types should be here.
 */
class InfoOut;
class InfoOut_buffer : public std::stringbuf {
  friend class InfoOut;

 protected:
  static const MInt m_fileFormatVersion = 1; //!< File format version (increase this by one every time you make
                                             //!< changes that could affect postprocessing tools)
  MBool m_rootOnly;                          //!< Stores whether only the root domain writes a log file
  MInt m_domainId;                           //!< Contains the MPI rank (= domain id) of this process
  MInt m_noDomains;                          //!< Contains the MPI rank count (= number of domains)
  MInt m_minFlushSize;                       //!< Minimum length of the internal buffer before flushing
  MString m_prefixMessage;                   //!< Stores the prefix that is prepended to each output
  MString m_suffixMessage;                   //!< Stores the suffix that is apended to each output
  std::ostringstream m_tmpBuffer;            //!< Temporary buffer to hold string until flushing
  MString m_projectName;                     //!< Name of the current Project

  std::vector<std::pair<MString, MString>> m_prefixAttributes;

  virtual MString encodeXml(const std::string& str);
  virtual MString getXmlHeader();
  virtual MString getXmlFooter();
  virtual void createPrefixMessage();
  virtual void createSuffixMessage();
  virtual void flushBuffer() = 0;

 public:
  InfoOut_buffer();
  virtual MBool setRootOnly(MBool rootOnly = true);
  virtual MInt setMinFlushSize(MInt minFlushSize);
};

/**
 * \brief Base class for all InfoOut<xyz> classes.
 * \author Michael Schlottke
 * \date June 2012
 * \details This class is used to hold stream/buffer-independent methods. All InfoOut<xyz> subclasses inherit
 *          from this class. The auxiliary classes InfoOut_<xyz> (especially the buffers), however, should NOT
 *          have this class as their parent.
 */
class InfoOut : public std::ostream {
  friend class InfoOut_buffer;

 protected:
  InfoOut_buffer* m_buffer = nullptr;

 public:
// Disable clang warning for ostream(m_buffer), although this is probably a compiler error
#if defined(MAIA_CLANG_COMPILER)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wuninitialized"
#endif
  InfoOut() : std::ostream(m_buffer){};
#if defined(MAIA_CLANG_COMPILER)
#pragma clang diagnostic pop
#endif
  virtual MBool setRootOnly(MBool rootOnly = true) = 0;
  MInt addAttribute(std::pair<MString, MString>);
  void eraseAttribute(MInt);
  void modifyAttribute(MInt, std::pair<MString, MString>);
};

/**
 * \brief Customized buffer to facilitate MPI I/O usage for a single file for all domains within an MPI communicator.
 * \author Michael Schlottke
 * \date June 2012
 * \details This class can be used as a regular string buffer, as it inherits from stringbuf. On flushing the
 *          buffer, the contents of the buffer are written to a file using MPI I/O. The entire MPI communication
 *          is hidden from the user, so that the underlying algorithms can be changed/optimized without affecting
 *          current implementations.
 *
 *          Internally, an XML file is created to store all messages with each processor's MPI rank so that each
 *          message can be attributed to the origin. The XML format ensures that fast and easy-to-develop post
 *          processing tools may be used (e.g. extractdomainlog.py). Additionally, meta information about the
 *          total number of domains as well as the creation/closing date of the file is saved.
 */
class InfoOut_mpiFileBuffer : public InfoOut_buffer {
 private:
  static const MInt m_maxStringLength = 8192; //!< Maximum string length (including formatting, default: 4096)
  MInt m_maxMessageLength;                    //!< Maximum message length (excluding formatting)
  MInt m_mpiWriteBufferSize;                  //!< Size of the MPI write buffer
  MChar* m_mpiWriteBuffer;                    //!< MPI write buffer
  MBool m_isOpen;                             //!< Stores whether the MPI file was already opened
  MString m_filename;                         //!< Filename on disk
  MPI_Comm m_mpiComm;                         //!< MPI communicator group
  MPI_File m_mpiFileHandle;                   //!< MPI file handle
  MPI_Request m_mpiRequest;                   //!< MPI request object (nonblocking I/O)


  MInt setMpiWriteBuffer(MInt newBufferSize);

 protected:
  virtual MInt sync();
  virtual void flushBuffer();

 public:
  InfoOut_mpiFileBuffer();
  InfoOut_mpiFileBuffer(const MString& filename, const MString& projectName, MPI_Comm mpiComm = globalMaiaCommWorld());
  ~InfoOut_mpiFileBuffer();
  void open(const MString& filename, const MString& projectName, MPI_Comm mpiComm = globalMaiaCommWorld());
  void close(MBool forceClose = false);
  MInt setMinFlushSize(MInt minFlushSize);
};

/**
 * \brief Customized buffer to facilitate of a regular physical file for each processor within an MPI communicator.
 * \author Michael Schlottke
 * \date June 2012
 * \details This class can be used as a regular string buffer, as it inherits from stringbuf. On flushing the
 *          buffer, the contents of the buffer are written to a file using an ofstream. This is mainly for cases
 *          where logging speed is crucial, as the implementation is very lightweight and since for each process
 *          an individual file is maintained. There is an option to use this buffer but to create only one file
 *          for the MPI root domain (see #m_rootOnlyHardwired).
 */
class InfoOut_simpleFileBuffer : public InfoOut_buffer {
 private:
  MBool m_isOpen;            //!< Stores whether the file(s) were already opened
  MBool m_rootOnlyHardwired; //!< If true, only domain 0 opens and uses a file
  MString m_filename;        //!< Filename on disk
  std::ofstream m_file;      //!< File stream tied to physical file on disk
  MPI_Comm m_mpiComm;        //!< MPI communicator group

 protected:
  virtual MInt sync();
  virtual void flushBuffer();

 public:
  InfoOut_simpleFileBuffer();
  InfoOut_simpleFileBuffer(const MString& filename, const MString& projectName, MPI_Comm mpiComm = globalMaiaCommWorld(),
                           MBool rootOnlyHardwired = false);
  ~InfoOut_simpleFileBuffer();
  void open(const MString& filename, const MString& projectName, MPI_Comm mpiComm = globalMaiaCommWorld(),
            MBool rootOnlyHardwired = false);
  void close(MBool forceClose = false);
};


/**
 * \brief Customized string buffer to prepend cout/cerr with the domain id.
 * \author Michael Schlottke
 * \date AJune 2012
 * \details This class takes an exisiting output stream such as cout or cerr and adds an additonal layer of
 *          buffering. For cerr this effectively means that cerr is now thread-safe, i.e. multiple processors
 *          writing to cerr at once should not result in mixed up strings anymore. Furthermore, each message is
 *          prepended with a prefix consisting of the rank of the specified MPI communicator. This behavior is
 *          optional and may be changed at runtime by subsequent calls to initialize().
 */
class InfoOut_streamBuffer : public InfoOut_buffer {
 private:
  MBool m_isInitialized;  //!< Stores whether a stream was already associated with this buffer
  MBool m_printDomainId;  //!< Stores whether the domain id should be prepended to each output
  std::ostream* m_output; //!< Stores the output stream (usually cout or cerr) that will be modified
  MPI_Comm m_mpiComm;     //!< MPI communicator group
  MBool m_isDisabled;     //!< Stores wether output in streamBuffer should be written

 protected:
  virtual MInt sync();
  virtual void flushBuffer();
  virtual MString addPrefix(const MString& str);

 public:
  InfoOut_streamBuffer();
  InfoOut_streamBuffer(std::ostream* os, MPI_Comm mpiComm = globalMaiaCommWorld(), MBool printDomainId = true);
  void initialize(std::ostream* os, MPI_Comm mpiComm = globalMaiaCommWorld(), MBool printDomainId = true);
};

/**
 * \brief Class to create a create an output stream for a writable file, using either MPI I/O or a physical file.
 * \author Michael Schlottke
 * \date June 2012
 * \details This class can be used to open a file on all processors (alternatively: only on a specified MPI
 *          communicator) and write to it using the normal C++ stream syntax (i.e. just like cout or cerr).
 *          Internally, it uses a InfoOut_mpiFileBuffer object as the internal buffer, and can thus hide all
 *          MPI-related commands from the user.
 *
 *          Alternatively, it is also possible to use regular physical files for each processor, and to write
 *          directly to it using a regular ofstream. This mode also allows for the setting that only process 0
 *          within an MPI communicator opens a file to write to.
 */
class InfoOutFile : public InfoOut {
 private:
  MInt m_fileType; //!< File type that is opened
  MBool m_isOpen;  //!< Stores whether a file was already opened or not

 public:
  InfoOutFile();
  InfoOutFile(const MString& filename, const MString& projectName, MInt fileType = 0, MPI_Comm mpiComm = globalMaiaCommWorld(),
              MBool rootOnlyHardwired = false);
  ~InfoOutFile();
  void open(const MString& filename, const MString& projectName, MInt fileType = 0, MPI_Comm mpiComm = globalMaiaCommWorld(),
            MBool rootOnlyHardwired = false);
  void close(MBool forceClose = false);
  MBool setRootOnly(MBool rootOnly = true);
  MInt setMinFlushSize(MInt minFlushSize);
};

/**
 * \brief Class to create an output stream that writes to cout or cerr but prepends each line with the MPI rank.
 * \author Michael Schlottke
 * \date June 2012
 * \details This class can be used to write to an already existing output stream. The usage and behavior is just
 *          like with existing C++ output streams such as cout/cerr, except that it offers the possibility to
 *          prepend each line with information about the current MPI domain id. In future implementations it might
 *          be worth considering to add even more information to each message.
 */
class InfoOutStream : public InfoOut {
 private:
  MBool m_isInitialized; //!< Stores whether a stream was already opened or not

 public:
  InfoOutStream();
  InfoOutStream(std::ostream* os, MPI_Comm mpiComm = globalMaiaCommWorld(), MBool printDomainId = true);
  ~InfoOutStream();
  void initialize(std::ostream* os, MPI_Comm mpiComm = globalMaiaCommWorld(), MBool printDomainId = true);
  MBool setRootOnly(MBool rootOnly = true);
};

#endif /* ifndef INFOOUT_H */
