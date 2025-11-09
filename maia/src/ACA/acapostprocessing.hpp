#include "acapostprocessing.h"

#include <vector>
#include "INCLUDE/maiatypes.h"

/** \brief ACA post processing class for overall sound pressure calculation
 *
 */
template <MInt nDim>
class AcaPostProcessingOASPL : public AcaPostProcessing {
 private:
  MString m_fileNameIoParallel;
  MFloat m_pRefSq = -1.0;
  const MFloat m_gamma = 1.4; /// isentropic exponent
  std::vector<MFloat> m_oaspl;

 public:
  using AcaPostProcessing::AcaPostProcessing;
  void init_() override;
  void calc(const MInt observerId, const MFloat* const data, const MFloat* const dataComplex = nullptr) override;
  void finish() override;
};

template <MInt nDim>
void AcaPostProcessingOASPL<nDim>::init_() {
  m_fileNameIoParallel = m_outPath + "OASPL" + ParallelIo::fileExt();
  m_oaspl.resize(m_noObservers, 0.0);

  constexpr MFloat pRefDim = 2e-5; // [Pa]
  MFloat pInftyDim = 101325;       // [Pa] = (rho_infty * a_infty^2 ) / gamma
  pInftyDim = Context::getBasicProperty<MFloat>("pInftyDimensional", AT_, &pInftyDim);
  m_pRefSq = POW2(pRefDim / (m_gamma * pInftyDim)); // [-]
}

template <MInt nDim>
void AcaPostProcessingOASPL<nDim>::calc(const MInt observerId, const MFloat* const data,
                                        const MFloat* const /*dataComplex*/) {
  // Calculation of overall sound pressure level using p_rms (could also be computed by the sum of
  // spl over all frequencies)
  MFloat p_rms = 0.0;
  for(MInt t = 0; t < m_noSamples; t++) {
    p_rms += POW2(data[t]);
  }
  p_rms = p_rms / m_noSamples;
  // Calculate overall sound pressure level for this observer
  m_oaspl[observerId] = 10.0 * log10(p_rms / m_pRefSq);
}

template <MInt nDim>
void AcaPostProcessingOASPL<nDim>::finish() {
  using namespace maia::parallel_io;
  ParallelIo file(m_fileNameIoParallel, PIO_REPLACE, m_mpiComm);
  ParallelIo::size_type dimSizesCoordinates[] = {m_noGlobalObservers, nDim};
  file.defineArray(PIO_FLOAT, "coordinates", 2, &dimSizesCoordinates[0]);
  file.defineArray(PIO_FLOAT, "OASPL", m_noGlobalObservers);

  const MBool isRoot = (m_rank == 0);
  file.setOffset((isRoot) ? m_noGlobalObservers : 0, 0, 2);
  file.writeArray(m_coords, "coordinates");
  file.setOffset(m_noObservers, m_offsetObserver);
  file.writeArray(m_oaspl.data(), "OASPL");
}

/** \brief ACA post processing class for absoulte pressure calculation
 *
 */
template <MInt nDim>
class AcaPostProcessingPABS : public AcaPostProcessing {
 private:
  MString m_fileNameIoParallel;
  std::vector<MFloat> m_pabs;

 public:
  using AcaPostProcessing::AcaPostProcessing;
  void init_() override;
  void calc(const MInt observerId, const MFloat* const data, const MFloat* const dataComplex = nullptr) override;
  void finish() override;
};

template <MInt nDim>
void AcaPostProcessingPABS<nDim>::init_() {
  m_fileNameIoParallel = m_outPath + "pressABS" + ParallelIo::fileExt();
  m_pabs.resize(m_noObservers * m_noSamples, 0.0);
}

template <MInt nDim>
void AcaPostProcessingPABS<nDim>::calc(const MInt observerId, const MFloat* const /*data*/,
                                       const MFloat* const dataComplex) {
  const MInt offset = observerId * m_noSamples;
  for(MInt t = 1; t < m_noSamples; t++) { // Skip frequency 0
    const MFloat real = dataComplex[t * 2 + 0];
    const MFloat imag = dataComplex[t * 2 + 1];
    // product of p * p_complex_conjugate = p_re_sq + p_im_sq = [abs(p)]^2
    const MFloat p_magnSq = POW2(real) + POW2(imag);
    // Absolute pressure
    m_pabs[offset + t] = 2 * sqrt(p_magnSq);
  }
}

template <MInt nDim>
void AcaPostProcessingPABS<nDim>::finish() {
  using namespace maia::parallel_io;
  ParallelIo file(m_fileNameIoParallel, PIO_REPLACE, m_mpiComm);
  ParallelIo::size_type dimSizesCoordinates[] = {m_noGlobalObservers, nDim};
  ParallelIo::size_type dimSizes[] = {m_noGlobalObservers, m_noSamples};
  file.defineArray(PIO_FLOAT, "coordinates", 2, &dimSizesCoordinates[0]);
  file.defineArray(PIO_FLOAT, "pressABS", 2, &dimSizes[0]);

  const MBool isRoot = (m_rank == 0);
  file.setOffset((isRoot) ? m_noGlobalObservers : 0, 0, 2);
  file.writeArray(m_coords, "coordinates");
  ScratchSpace<MFloat> pressABS(m_noObservers, m_noSamples, AT_, "pressABS");
  file.setOffset(m_noObservers, m_offsetObserver, 2);
  for(MInt obs = 0; obs < m_noObservers; obs++) {
    const MInt offset = obs * m_noSamples;
    for(MInt t = 0; t < m_noSamples; t++) {
      pressABS(obs, t) = m_pabs[offset + t];
    }
  }
  file.writeArray(pressABS.data(), "pressABS");
}

/** \brief ACA post processing class for calculation of root mean square pressure
 *
 */
template <MInt nDim>
class AcaPostProcessingRMS : public AcaPostProcessing {
 private:
  MString m_fileNameIoParallel;
  MString m_fileNameIoDat;
  std::vector<MFloat> m_prms;

 public:
  using AcaPostProcessing::AcaPostProcessing;
  void init_() override;
  void calc(const MInt observerId, const MFloat* const data, const MFloat* const dataComplex = nullptr) override;
  void finish() override;
};

template <MInt nDim>
void AcaPostProcessingRMS<nDim>::init_() {
  m_fileNameIoParallel = m_outPath + "RMSpressure" + ParallelIo::fileExt();
  m_fileNameIoDat = m_outPath + "RMS.dat";
  m_prms.resize(m_noGlobalObservers, 0.0);
}

template <MInt nDim>
void AcaPostProcessingRMS<nDim>::calc(const MInt observerId, const MFloat* const data,
                                      const MFloat* const /*dataComplex*/) {
  // Calculation of root mean square pressure values
  m_prms[m_offsetObserver + observerId] = 0;
  for(MInt t = 0; t < m_noSamples; t++) {
    m_prms[m_offsetObserver + observerId] += POW2(data[t]);
  }
  m_prms[m_offsetObserver + observerId] = sqrt(m_prms[m_offsetObserver + observerId] / m_noSamples);
}

template <MInt nDim>
void AcaPostProcessingRMS<nDim>::finish() {
  constexpr MInt rootId = 0;
  using namespace maia::parallel_io;
  ParallelIo file(m_fileNameIoParallel, PIO_REPLACE, m_mpiComm);
  ParallelIo::size_type dimSizesCoordinates[] = {m_noGlobalObservers, nDim};
  file.defineArray(PIO_FLOAT, "coordinates", 2, &dimSizesCoordinates[0]);
  file.defineArray(PIO_FLOAT, "rms_pressure", m_noGlobalObservers);

  const MBool isRoot = (m_rank == rootId);
  file.setOffset((isRoot) ? m_noGlobalObservers : 0, 0, 2);
  file.writeArray(m_coords, "coordinates");
  if(isRoot) {
    MPI_Reduce(MPI_IN_PLACE, m_prms.data(), m_noGlobalObservers, maia::type_traits<MFloat>::mpiType(), MPI_SUM, rootId,
               m_mpiComm, AT_, "MPI_IN_PLACE", "m_prms");
    // TODO labels:ACA,totest Check distl, theta, and alpha
    std::ofstream outfile;
    outfile.open(m_fileNameIoDat);
    if constexpr(nDim == 2) {
      outfile << "# obs:1 x:2 y:3 theta:4 p_rms:5" << std::endl;
      for(MInt obs = 0; obs < m_noGlobalObservers; obs++) {
        const MFloat* const obsCoord = &m_coords[nDim * obs];
        const MFloat theta = atan2(obsCoord[1], obsCoord[0]);
        outfile << obs << " " << obsCoord[0] << " " << obsCoord[1] << " " << theta << " " << m_prms[obs] << std::endl;
      }
    } else if constexpr(nDim == 3) {
      outfile << "# obs:1 x:2 y:3 z:4 theta:5 alpha:6 p_rms:7" << std::endl;
      for(MInt obs = 0; obs < m_noGlobalObservers; obs++) {
        const MFloat* const obsCoord = &m_coords[nDim * obs];
        const MFloat dist = sqrt(POW2(obsCoord[0]) + POW2(obsCoord[1]) + POW2(obsCoord[2]));
        const MFloat theta = atan2(obsCoord[2], obsCoord[0]);
        const MFloat alpha = asin(obsCoord[1] / dist);
        outfile << obs << " " << obsCoord[0] << " " << obsCoord[1] << " " << obsCoord[2] << " " << theta << " " << alpha
                << " " << m_prms[obs] << std::endl;
      }
    }
    outfile.close();
    file.setOffset(m_noGlobalObservers, 0);
  } else {
    MPI_Reduce(m_prms.data(), nullptr, m_noGlobalObservers, maia::type_traits<MFloat>::mpiType(), MPI_SUM, rootId,
               m_mpiComm, AT_, "MPI_IN_PLACE", "m_prms");
    file.setOffset(0, 0);
  }
  file.writeArray(m_prms.data(), "rms_pressure");
}

/** \brief ACA post processing class for sound pressure level calculation
 *
 */
template <MInt nDim>
class AcaPostProcessingSPL : public AcaPostProcessing {
 private:
  MString m_fileNameIoParallel;
  MFloat m_pRefSq = -1.0;
  const MFloat m_gamma = 1.4; /// isentropic exponent
  std::vector<MFloat> m_spl;

 public:
  using AcaPostProcessing::AcaPostProcessing;
  void init_() override;
  void calc(const MInt observerId, const MFloat* const data, const MFloat* const dataComplex = nullptr) override;
  void finish() override;
};

template <MInt nDim>
void AcaPostProcessingSPL<nDim>::init_() {
  m_fileNameIoParallel = m_outPath + "SPL" + ParallelIo::fileExt();
  // reference pressure: Attention: The default value needs to be made non-dimensional
  /*! \property
    \page propertyPageACA ACA
    \section pInftyDimensional
    <code>MFloat pInftyDim</code>\n
    default = <code>101325</code>\n\n
    Used to calculate non-dimensional reference pressure. Is given in [Pa] by
    the expression (rho_infty * a_infty^2)/gamma. This is needed, as p_ref is
    defined to be 2e-5 Pa for air.
    \n\n
    Keywords: <i>ACA, post</i>
  */
  constexpr MFloat pRefDim = 2e-5; // [Pa]
  MFloat pInftyDim = 101325;       // [Pa] = (rho_infty * a_infty^2 ) / gamma
  pInftyDim = Context::getBasicProperty<MFloat>("pInftyDimensional", AT_, &pInftyDim);
  m_pRefSq = POW2(pRefDim / (m_gamma * pInftyDim)); // [-]

  m_spl.resize(m_noObservers * m_noSamples, 0.0);
}

template <MInt nDim>
void AcaPostProcessingSPL<nDim>::calc(const MInt observerId, const MFloat* const /*data*/,
                                      const MFloat* const dataComplex) {
  // Calculation of sound pressure level (pressure at each frequency)
  const MInt offset = observerId * m_noSamples;
  for(MInt t = 1; t < m_noSamples; t++) { // Skip frequency 0
    const MFloat real = dataComplex[t * 2 + 0];
    const MFloat imag = dataComplex[t * 2 + 1];
    // p times p complex conjugate = p_re_sq + p_imag_sq
    const MFloat pSq = POW2(real) + POW2(imag);
    // Sound pressure for each observer point for each frequency
    // Ref.:  Mendez et al. (2013): https://doi.org/10.1260%2F1475-472X.12.1-2.1
    //        or Beranek and Ver (1992); Delfs, Basics of Aeroacoustics (2016); ...
    // Remind: 20 * log(a) = 10 * log(a^2), and pSq is only halved value because of FFT
    m_spl[offset + t] = 10 * log10(2 * pSq / m_pRefSq);
  }
}

template <MInt nDim>
void AcaPostProcessingSPL<nDim>::finish() {
  using namespace maia::parallel_io;
  ParallelIo file(m_fileNameIoParallel, PIO_REPLACE, m_mpiComm);
  ParallelIo::size_type dimSizesCoordinates[] = {m_noGlobalObservers, nDim};
  ParallelIo::size_type dimSizesSpl[] = {m_noGlobalObservers, m_noSamples};
  file.defineArray(PIO_FLOAT, "coordinates", 2, &dimSizesCoordinates[0]);
  file.defineArray(PIO_FLOAT, "frequency", m_noSamples);
  file.defineArray(PIO_FLOAT, "SPL", 2, &dimSizesSpl[0]);

  const MBool isRoot = (m_rank == 0);
  file.setOffset((isRoot) ? m_noGlobalObservers : 0, 0, 2);
  file.writeArray(m_coords, "coordinates");
  file.setOffset((isRoot) ? m_noSamples : 0, 0);
  file.writeArray(m_frequencies, "frequency");
  file.setOffset(m_noObservers, m_offsetObserver, 2);
  file.writeArray(m_spl.data(), "SPL");
}
