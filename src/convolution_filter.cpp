#include "../include/convolution_filter.h"

#include <assert.h>
#include <string.h>

#if defined(__APPLE__)
#include <Accelerate/Accelerate.h>
#endif

ConvolutionFilter::ConvolutionFilter() {
    m_shiftRegister = nullptr;
    m_impulseResponse = nullptr;

    m_shiftOffset = 0;
    m_sampleCount = 0;
}

ConvolutionFilter::~ConvolutionFilter() {
    assert(m_shiftRegister == nullptr);
    assert(m_impulseResponse == nullptr);
}

void ConvolutionFilter::initialize(int samples) {
    m_sampleCount = samples;
    m_shiftOffset = 0;
    m_shiftRegister = new float[samples];
    m_impulseResponse = new float[samples];

    memset(m_shiftRegister, 0, sizeof(float) * samples);
    memset(m_impulseResponse, 0, sizeof(float) * samples);
}

void ConvolutionFilter::destroy() {
    delete[] m_shiftRegister;
    delete[] m_impulseResponse;
    m_shiftRegister = nullptr;
    m_impulseResponse = nullptr;
#if defined(__APPLE__)
    if (m_fftSetup != nullptr) vDSP_destroy_fftsetup(static_cast<FFTSetup>(m_fftSetup));
    m_fftSetup = nullptr;
    delete[] m_hRe;
    delete[] m_hIm;
    delete[] m_xRe;
    delete[] m_xIm;
    m_hRe = m_hIm = m_xRe = m_xIm = nullptr;
    m_fftTried = m_fftReady = false;
#endif
}

#if defined(__APPLE__)
// nfs-mac: called on the first sample (the synthesizer fills the impulse response after initialize()).
bool ConvolutionFilter::prepareFft() {
    m_fftTried = true;
    if (m_sampleCount < 2 * kBlock) return false;
    constexpr int log2n = 8;  // FFT size 2 x kBlock = 256
    static_assert((1 << log2n) == 2 * kBlock, "FFT size");
    FFTSetup setup = vDSP_create_fftsetup(log2n, FFT_RADIX2);
    if (setup == nullptr) return false;
    m_fftSetup = setup;
    m_partitions = (m_sampleCount + kBlock - 1) / kBlock;
    const size_t n = size_t(m_partitions) * kBlock;
    m_hRe = new float[n];
    m_hIm = new float[n];
    m_xRe = new float[n]();
    m_xIm = new float[n]();
    for (int p = 0; p < m_partitions; ++p) {
        memset(m_tmp, 0, sizeof m_tmp);
        const int left = m_sampleCount - p * kBlock;
        const int count = left < kBlock ? left : kBlock;
        memcpy(m_tmp, m_impulseResponse + p * kBlock, sizeof(float) * size_t(count));
        DSPSplitComplex h = {m_hRe + p * kBlock, m_hIm + p * kBlock};
        vDSP_ctoz(reinterpret_cast<const DSPComplex *>(m_tmp), 2, &h, 1, kBlock);
        vDSP_fft_zrip(setup, &h, 1, log2n, FFT_FORWARD);
    }
    memset(m_in, 0, sizeof m_in);
    memset(m_out, 0, sizeof m_out);
    m_fdlPos = 0;
    m_pos = 0;
    return true;
}

// nfs-mac: one block through every partition (overlap-save: the last kBlock samples of the inverse are valid).
void ConvolutionFilter::processBlock() {
    constexpr int log2n = 8;
    FFTSetup setup = static_cast<FFTSetup>(m_fftSetup);
    DSPSplitComplex x = {m_xRe + m_fdlPos * kBlock, m_xIm + m_fdlPos * kBlock};
    vDSP_ctoz(reinterpret_cast<const DSPComplex *>(m_in), 2, &x, 1, kBlock);
    vDSP_fft_zrip(setup, &x, 1, log2n, FFT_FORWARD);

    memset(m_accRe, 0, sizeof m_accRe);
    memset(m_accIm, 0, sizeof m_accIm);
    float dc = 0, nyquist = 0;  // packed: element 0 holds DC (real part) and Nyquist (imaginary part), both real
    DSPSplitComplex acc = {m_accRe + 1, m_accIm + 1};
    for (int p = 0; p < m_partitions; ++p) {
        const int q = (m_fdlPos - p + m_partitions) % m_partitions;  // the input block p blocks ago
        float *xr = m_xRe + q * kBlock, *xi = m_xIm + q * kBlock;
        float *hr = m_hRe + p * kBlock, *hi = m_hIm + p * kBlock;
        dc += xr[0] * hr[0];
        nyquist += xi[0] * hi[0];
        DSPSplitComplex xs = {xr + 1, xi + 1}, hs = {hr + 1, hi + 1};
        vDSP_zvma(&xs, 1, &hs, 1, &acc, 1, &acc, 1, kBlock - 1);
    }
    m_accRe[0] = dc;
    m_accIm[0] = nyquist;
    DSPSplitComplex y = {m_accRe, m_accIm};
    vDSP_fft_zrip(setup, &y, 1, log2n, FFT_INVERSE);
    vDSP_ztoc(&y, 1, reinterpret_cast<DSPComplex *>(m_tmp), 2, kBlock);
    // vDSP's forward real FFT is 2x the DFT (input and IR: 4x) and its inverse N x: 1 / (4 N), N = 2 kBlock.
    const float scale = 1.0f / (8.0f * kBlock);
    vDSP_vsmul(m_tmp + kBlock, 1, &scale, m_out, 1, kBlock);

    m_fdlPos = (m_fdlPos + 1) % m_partitions;
    memcpy(m_in, m_in + kBlock, sizeof(float) * kBlock);
}
#endif

float ConvolutionFilter::f(float sample) {
#if defined(__APPLE__)
    if (!m_fftTried) m_fftReady = prepareFft();
    if (m_fftReady) {  // nfs-mac: the FFT form (kBlock samples of latency)
        m_in[kBlock + m_pos] = sample;
        const float result = m_out[m_pos];
        if (++m_pos == kBlock) {
            processBlock();
            m_pos = 0;
        }
        return result;
    }
#endif
    m_shiftRegister[m_shiftOffset] = sample;

    float result = 0;
    for (int i = 0; i < m_sampleCount - m_shiftOffset; ++i) {
        result += m_impulseResponse[i] * m_shiftRegister[i + m_shiftOffset];
    }

    for (int i = m_sampleCount - m_shiftOffset; i < m_sampleCount; ++i) {
        result += m_impulseResponse[i] * m_shiftRegister[i - (m_sampleCount - m_shiftOffset)];
    }

    m_shiftOffset = (m_shiftOffset - 1 + m_sampleCount) % m_sampleCount;

    return result;
}
