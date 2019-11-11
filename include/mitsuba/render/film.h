#pragma once

#include <mitsuba/mitsuba.h>
#include <mitsuba/core/logger.h>
#include <mitsuba/core/object.h>
#include <mitsuba/core/rfilter.h>
#include <mitsuba/core/vector.h>
#include <mitsuba/render/sampler.h>
#include <mitsuba/render/fwd.h>

NAMESPACE_BEGIN(mitsuba)

/** \brief Abstract film base class - used to store samples
 * generated by \ref Integrator implementations.
 *
 * To avoid lock-related bottlenecks when rendering with many cores,
 * rendering threads first store results in an "image block", which
 * is then committed to the film using the \ref put() method.
 */
template <typename Float, typename Spectrum>
class MTS_EXPORT_RENDER Film : public Object {
public:
    MTS_DECLARE_CLASS_VARIANT(Film, Object, "film")
    MTS_IMPORT_TYPES();
    using ImageBlock           = typename RenderAliases::ImageBlock;
    using ReconstructionFilter = ReconstructionFilter<Float, Spectrum>;

    /// Clear the film
    virtual void clear() = 0;

    /// Merge an image block into the film
    virtual void put(const ImageBlock *block) = 0;

    /// Overwrite the film with the given bitmap.
    /// The size of the given bitmap has to match to current size.
    virtual void set_bitmap(const Bitmap *bitmap) = 0;

    /// Accumulate a bitmap on top of the radiance values stored in the film.
    /// The size of the given bitmap has to match to current size.
    virtual void add_bitmap(const Bitmap *bitmap, ScalarFloat multiplier = 1.0f) = 0;

    /// Set the target filename (with or without extension)
    virtual void set_destination_file(const fs::path &filename) = 0;

    /// Develop the film and write the result to the previously specified filename
    virtual void develop() = 0;

    /**
     * \brief Develop the contents of a subregion of the film and store
     * it inside the given bitmap
     *
     * This may fail when the film does not have an explicit representation
     * of the bitmap in question (e.g. when it is writing to a tiled EXR image)
     *
     * \return \c true upon success
     */
    virtual bool develop(
        const ScalarPoint2i  &offset,
        const ScalarVector2i &size,
        const ScalarPoint2i  &target_offset,
        Bitmap *target) const = 0;

    /// Does the destination file already exist?
    virtual bool destination_exists(const fs::path &basename) const = 0;

    /**
     * Should regions slightly outside the image plane be sampled to improve
     * the quality of the reconstruction at the edges? This only makes
     * sense when reconstruction filters other than the box filter are used.
     */
    bool has_high_quality_edges() const { return m_high_quality_edges; }

    // =============================================================
    //! @{ \name Accessor functions
    // =============================================================

    /// Returns a pointer to the underlying image storage, or nullptr
    /// if there is none.
    virtual Bitmap *bitmap() { return nullptr; }

    /// Ignoring the crop window, return the resolution of the underlying sensor
    const ScalarVector2i &size() const { return m_size; }

    /// Return the size of the crop window
    const ScalarVector2i &crop_size() const { return m_crop_size; }

    /// Return the offset of the crop window
    const ScalarPoint2i &crop_offset() const { return m_crop_offset; }

    /**
     * Set the size and offset of the crop window.
     * Values stored in the film may be reset when crop size changes.
     */
    virtual void set_crop_window(const ScalarVector2i &crop_size, const ScalarPoint2i &crop_offset);

    /// Return the image reconstruction filter (const version)
    const ReconstructionFilter *reconstruction_filter() const {
        return m_filter.get();
    }

    //! @}
    // =============================================================

    virtual std::string to_string() const override {
        std::ostringstream oss;
        oss << "Film[" << std::endl
            << "  size = "        << m_size        << "," << std::endl
            << "  crop_size = "   << m_crop_size   << "," << std::endl
            << "  crop_offset = " << m_crop_offset << "," << std::endl
            << "  high_quality_edges = " << m_high_quality_edges << "," << std::endl
            << "  m_filter = " << m_filter << std::endl
            << "]";
        return oss.str();
    }

protected:
    /// Create a film
    Film(const Properties &props);

    /// Virtual destructor
    virtual ~Film();

    /// Throws if the crop window specification is invalid.
    void check_valid_crop_window() const;

protected:
    ScalarVector2i m_size, m_crop_size;
    ScalarPoint2i  m_crop_offset;
    bool m_high_quality_edges;
    ref<ReconstructionFilter> m_filter;
};

NAMESPACE_END(mitsuba)
