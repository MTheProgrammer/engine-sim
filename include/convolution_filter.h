#ifndef ATG_ENGINE_SIM_CONVOLUTION_FILTER_H
#define ATG_ENGINE_SIM_CONVOLUTION_FILTER_H

#include "filter.h"

class ConvolutionFilter : public Filter {
    public:
        ConvolutionFilter();
        virtual ~ConvolutionFilter();

        void initialize(int samples);
        virtual float f(float sample) override;
        virtual void destroy();

        int getSampleCount() const { return m_sampleCount; }
        float *getImpulseResponse() { return m_impulseResponse; }

    protected:
        float *m_shiftRegister;
        int m_shiftOffset;

        float *m_impulseResponse;
        int m_sampleCount;

#if defined(__APPLE__)
        // nfs-mac: uniformly partitioned FFT convolution (overlap-save, Accelerate vDSP), blocks of kBlock samples.
        // The direct form above cost one multiply-add per tap per sample (~29 % of a Viper V10's CPU); this is about
        // 50x less. The output is kBlock samples (2.9 ms) later. Impulse responses under 2 x kBlock taps keep the
        // direct form.
        static constexpr int kBlock = 128;
        bool prepareFft();
        void processBlock();
        bool m_fftTried = false;
        bool m_fftReady = false;
        void *m_fftSetup = nullptr;  // FFTSetup
        int m_partitions = 0;
        float *m_hRe = nullptr, *m_hIm = nullptr;  // the IR's partitions' spectra, m_partitions x kBlock (packed)
        float *m_xRe = nullptr, *m_xIm = nullptr;  // the input blocks' spectra (frequency-domain delay line)
        int m_fdlPos = 0;
        int m_pos = 0;
        float m_in[2 * kBlock];  // the previous block, then the block being filled
        float m_out[kBlock];     // the last processed block's output, read out while the next block fills
        float m_accRe[kBlock], m_accIm[kBlock], m_tmp[2 * kBlock];
#endif
};

#endif /* ATG_ENGINE_SIM_CONVOLUTION_FILTER_H */
