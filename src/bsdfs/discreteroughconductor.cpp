#include <mitsuba/core/string.h>
#include <mitsuba/core/fwd.h>
#include <mitsuba/core/plugin.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/ior.h>
#include <mitsuba/render/microfacet.h>
#include <mitsuba/render/texture.h>

#include <cmath>
#include <queue>

NAMESPACE_BEGIN(mitsuba)

/**!
.. _bsdf-roughconductor:

このプログラムはroughconductorのコピーを基に作っている。
Rough conductor material (:monosp:`roughconductor`)
---------------------------------------------------

.. pluginparameters::

 * - material
   - |string|
   - Name of the material preset, see :num:`conductor-ior-list`. (Default: none)

 * - eta, k
   - |spectrum| or |texture|
   - Real and imaginary components of the material's index of refraction. (Default: based on the value of :monosp:`material`)
   - |exposed|, |differentiable|, |discontinuous|

 * - specular_reflectance
   - |spectrum| or |texture|
   - Optional factor that can be used to modulate the specular reflection component.
     Note that for physical realism, this parameter should never be touched. (Default: 1.0)
   - |exposed|, |differentiable|

 * - distribution
   - |string|
   - Specifies the type of microfacet normal distribution used to model the surface roughness.

     - :monosp:`beckmann`: Physically-based distribution derived from Gaussian random surfaces.
       This is the default.
     - :monosp:`ggx`: The GGX :cite:`Walter07Microfacet` distribution (also known as Trowbridge-Reitz
       :cite:`Trowbridge19975Average` distribution) was designed to better approximate the long
       tails observed in measurements of ground surfaces, which are not modeled by the Beckmann
       distribution.

 * - alpha, alpha_u, alpha_v
   - |texture| or |float|
   - Specifies the roughness of the unresolved surface micro-geometry along the tangent and
     bitangent directions. When the Beckmann distribution is used, this parameter is equal to the
     **root mean square** (RMS) slope of the microfacets. :monosp:`alpha` is a convenience
     parameter to initialize both :monosp:`alpha_u` and :monosp:`alpha_v` to the same value. (Default: 0.1)
   - |exposed|, |differentiable|, |discontinuous|

 * - sample_visible
   - |bool|
   - Enables a sampling technique proposed by Heitz and D'Eon :cite:`Heitz1014Importance`, which
     focuses computation on the visible parts of the microfacet normal distribution, considerably
     reducing variance in some cases. (Default: |true|, i.e. use visible normal sampling)

This plugin implements a realistic microfacet scattering model for rendering
rough conducting materials, such as metals.

.. subfigstart::
.. subfigure:: ../../resources/data/docs/images/render/bsdf_roughconductor_copper.jpg
   :caption: Rough copper (Beckmann, :math:`\alpha=0.1`)
.. subfigure:: ../../resources/data/docs/images/render/bsdf_roughconductor_anisotropic_aluminium.jpg
   :caption: Vertically brushed aluminium (Anisotropic Beckmann, :math:`\alpha_u=0.05,\ \alpha_v=0.3`)
.. subfigure:: ../../resources/data/docs/images/render/bsdf_roughconductor_textured_carbon.jpg
   :caption: Carbon fiber using two inverted checkerboard textures for ``alpha_u`` and ``alpha_v``
.. subfigend::
    :label: fig-bsdf-roughconductor


Microfacet theory describes rough surfaces as an arrangement of unresolved
and ideally specular facets, whose normal directions are given by a
specially chosen *microfacet distribution*. By accounting for shadowing
and masking effects between these facets, it is possible to reproduce the
important off-specular reflections peaks observed in real-world measurements
of such materials.

This plugin is essentially the *roughened* equivalent of the (smooth) plugin
:ref:`conductor <bsdf-conductor>`. For very low values of :math:`\alpha`, the two will
be identical, though scenes using this plugin will take longer to render
due to the additional computational burden of tracking surface roughness.

The implementation is based on the paper *Microfacet Models
for Refraction through Rough Surfaces* by Walter et al.
:cite:`Walter07Microfacet` and it supports two different types of microfacet
distributions.

To facilitate the tedious task of specifying spectrally-varying index of
refraction information, this plugin can access a set of measured materials
for which visible-spectrum information was publicly available
(see the corresponding table in the :ref:`conductor <bsdf-conductor>` reference).

When no parameters are given, the plugin activates the default settings,
which describe a 100% reflective mirror with a medium amount of roughness modeled
using a Beckmann distribution.

To get an intuition about the effect of the surface roughness parameter
:math:`\alpha`, consider the following approximate classification: a value of
:math:`\alpha=0.001-0.01` corresponds to a material with slight imperfections
on an otherwise smooth surface finish, :math:`\alpha=0.1` is relatively rough,
and :math:`\alpha=0.3-0.7` is **extremely** rough (e.g. an etched or ground
finish). Values significantly above that are probably not too realistic.


The following XML snippet describes a material definition for brushed aluminium:

.. tabs::
    .. code-tab:: xml
        :name: lst-roughconductor-aluminium

        <bsdf type="roughconductor">
            <string name="material" value="Al"/>
            <string name="distribution" value="ggx"/>
            <float name="alpha_u" value="0.05"/>
            <float name="alpha_v" value="0.3"/>
        </bsdf>

    .. code-tab:: python

        'type': 'roughconductor',
        'material': 'Al',
        'distribution': 'ggx',
        'alpha_u': 0.05,
        'alpha_v': 0.3

Technical details
*****************

All microfacet distributions allow the specification of two distinct
roughness values along the tangent and bitangent directions. This can be
used to provide a material with a *brushed* appearance. The alignment
of the anisotropy will follow the UV parameterization of the underlying
mesh. This means that such an anisotropic material cannot be applied to
triangle meshes that are missing texture coordinates.

Since Mitsuba 0.5.1, this plugin uses a new importance sampling technique
contributed by Eric Heitz and Eugene D'Eon, which restricts the sampling
domain to the set of visible (unmasked) microfacet normals. The previous
approach of sampling all normals is still available and can be enabled
by setting :monosp:`sample_visible` to :monosp:`false`. However this will lead
to significantly slower convergence.

When using this plugin, you should ideally compile Mitsuba with support for
spectral rendering to get the most accurate results. While it also works
in RGB mode, the computations will be more approximate in nature.
Also note that this material is one-sided---that is, observed from the
back side, it will be completely black. If this is undesirable,
consider using the :ref:`twosided <bsdf-twosided>` BRDF adapter.

In *polarized* rendering modes, the material automatically switches to a polarized
implementation of the underlying Fresnel equations.

 */


// テスクチャ空間上の正方形を表すクラス
class MySquare {
    public:
        // 正方形の左上の座標
        float x1, y1;
        // 正方形の右下の座標
        float x2, y2;
        // 正方形の辺の長さ
        float width;
        // 正方形の面積
        float area;

        // コンストラクタ
        MySquare(float X1, float Y1, float X2, float Y2) : x1(X1), y1(Y1), x2(X2), y2(Y2) {
            width = std::abs(x1 - x2);
            area = width * width;
        }

        // 面積を定義したコンストラクタ
        MySquare(float X1, float Y1, float X2, float Y2, float Area) : x1(X1), y1(Y1), x2(X2), y2(Y2), area(Area) {
            width = std::abs(x1 - x2);
        }

        // 入力がない時のコンストラクタ
        MySquare() : x1(0.0f), y1(0.0f), x2(1.0f), y2(1.0f) {
            width = 1.0f;
            area = 1.0f;
        }

        
        // 正方形を４つの子に分割する関数
        MySquare* split() {
            MySquare result[4];

            // 左上の子
            result[0] = MySquare(x1, y1, x1 + width/2, y1 + width/2, area/4);
            // 右上の子
            result[1] = MySquare(x1 + width/2, y1, x2, y1 + width/2, area/4);
            // 左下の子
            result[2] = MySquare(x1, y1 + width/2, x1 + width/2, y2, area/4);
            // 右下の子
            result[3] = MySquare(x1 + width/2, y1 + width/2, x2, y2, area/4);

            return result;
        }
};
// ------------------------------------------------------------

// 方向領域の球面上にある三角形のクラス
class MySpheTri {
    public:
        // 各頂点の方向ベクトル
        dr::Array<float, 3> w1;
        dr::Array<float, 3> w2;
        dr::Array<float, 3> w3;
        // 三角形の面積
        float area;

        // コンストラクタ
        MySpheTri(dr::Array<float, 3> W1, dr::Array<float, 3> W2, dr::Array<float, 3> W3, float Area) : w1(W1), w2(W2), w3(W3), area(Area) {}

        // 入力がない時のコンストラクタ
        MySpheTri() {
            dr::Array<float, 3> tmp(0.0f, 0.0f, 0.0f);
            w1 = tmp;
            w2 = tmp;
            w3 = tmp;
            area = 0.0f; 
        }

        // 三角形の面積を求める関数（ハッシュテーブルから値を参照する）
        // float area()

        // 球面上にある三角形を４つの子に分割
        MySpheTri* split() {
            MySpheTri result[4];

            // 頂点w1を含む子
            result[0] = MySpheTri(w1, (w1 + w2)/2, (w1 + w3)/2, area/4);
            // 頂点w2を含む子
            result[1] = MySpheTri((w2 + w1)/2, w2, (w2 + w3)/2, area/4);
            // 中点のみで構成される子
            result[2] = MySpheTri((w2 + w3)/2, (w1 + w2)/2, (w3 + w1)/2, area/4);
            // 頂点w3を含む子
            result[3] = MySpheTri((w3 + w1)/2, (w3 + w2)/2, w3, area/4);

            return result;
        }

};
// ----------------------------------------------------------------------

// nodeを表すクラス
class MyNode {
    public:
        MySquare squ;
        MySpheTri tri;
        int num_particle;
        float vol;

        // 正方形と三角形が与えられた場合のコンストラクタ
        MyNode(MySquare Square, MySpheTri Triangle, int N_Particle) : squ(Square), tri(Triangle), num_particle(N_Particle) {
            vol = squ.area * tri.area;
        }

        // 入力がない時のコンストラクタ
        MyNode() {
            squ = MySquare();
            tri = MySpheTri();
            num_particle = 0;
            vol = 0.0;
        }

};
// -----------------------------------------------------------------------------------

// 平行四辺形の
class MyParallelogram {
    public:
        float x1, x2, x3, x4;
        float y1, y2, y3, y4;
        
        // コンストラクタ
        MyParallelogram(float X1, float Y1, float X2, float Y2, float X3, float Y3, float X4, float Y4) : 
        x2(0.0f), y2(0.0f), x3(0.0f), y3(0.0f), x4(0.0f), y4(0.0f) {
            // x座標が一番小さい頂点を(x1, y1)にする
            // 初期値は(x1, y1)=(X1, Y1)
            x1 = X1;
            y1 = Y1;
            float Vs[4][2] = {{X1, Y1}, {X2, Y2}, {X3, Y3}, {X4, Y4}};
            // x座標が一番小さい頂点のインデックス
            int idx_minX = 0;

            for(int i=1; i < 4; i++) {
                x1 = x1 <= Vs[i][0] ? x1 : Vs[i][0];
                y1 = x1 <= Vs[i][0] ? y1 : Vs[i][1];
                idx_minX = x1 <= Vs[i][0] ? idx_minX : i;
            }

            // x座標が１番小さい頂点を基準に、頂点を反時計回りにソートする
            for(int i = 0; i < 4; i++) {
                if(i = idx_minX) {continue;}
                
                else {

                }
            }
            

        }
        
        // 入力がない時のコンストラクタ
        MyParallelogram() : x1(0.0f), y1(0.0f), x2(0.0f), y2(0.0f), x3(0.0f), y3(0.0f), x4(0.0f), y4(0.0f) {}

        // 正方形との交差を判定（正方形内部に平行四辺形の頂点があればTrue）
        bool isIntersect(MySquare squ) {
            bool result = false;

            // 正方形が存在する範囲
            float minX = squ.x1 < squ.x2 ? squ.x1 : squ.x2;
            float maxX = squ.x1 > squ.x2 ? squ.x1 : squ.x2;
            float minY = squ.y1 < squ.y2 ? squ.y1 : squ.y2;
            float maxY = squ.y1 > squ.y2 ? squ.y1 : squ.y2;

            // 平行四辺形の各頂点
            float V[4][2] = {{x1, y1}, {x2, y2}, {x3, y3}, {x4, y4}};


            // 平行四辺形の頂点のうち正方形の内部に存在するものがあれば、交差するとしてtrueを返す
            for (int i = 0; i < 4; i++) {
                if((minX <= V[i][0]) && (V[i][0] <= maxX) && (minY <= V[i][1]) && (V[i][1] <= maxY)) {
                    result = true;
                    break;
                }
            }

            return result;
        }
};

template <typename Float, typename Spectrum>
class DiscreteRoughConductor final : public BSDF<Float, Spectrum> {
public:
    MI_IMPORT_BASE(BSDF, m_flags, m_components)
    MI_IMPORT_TYPES(Texture, MicrofacetDistribution)

    DiscreteRoughConductor(const Properties &props) : Base(props) {
        std::string material = props.string("material", "none");
        if (props.has_property("eta") || material == "none") {
            m_eta = props.texture<Texture>("eta", 0.f);
            m_k   = props.texture<Texture>("k",   1.f);
            if (material != "none")
                Throw("Should specify either (eta, k) or material, not both.");
        } else {
            std::tie(m_eta, m_k) = complex_ior_from_file<Spectrum, Texture>(props.string("material", "Cu"));
        }

        if (props.has_property("distribution")) {
            std::string distr = string::to_lower(props.string("distribution"));
            if (distr == "beckmann")
                m_type = MicrofacetType::Beckmann;
            else if (distr == "ggx")
                m_type = MicrofacetType::GGX;
            else
                Throw("Specified an invalid distribution \"%s\", must be "
                      "\"beckmann\" or \"ggx\"!", distr.c_str());
        } else {
            m_type = MicrofacetType::Beckmann;
        }

        m_sample_visible = props.get<bool>("sample_visible", true);

        if (props.has_property("alpha_u") || props.has_property("alpha_v")) {
            if (!props.has_property("alpha_u") || !props.has_property("alpha_v"))
                Throw("Microfacet model: both 'alpha_u' and 'alpha_v' must be specified.");
            if (props.has_property("alpha"))
                Throw("Microfacet model: please specify"
                      "either 'alpha' or 'alpha_u'/'alpha_v'.");
            m_alpha_u = props.texture<Texture>("alpha_u");
            m_alpha_v = props.texture<Texture>("alpha_v");
        } else {
            m_alpha_u = m_alpha_v = props.texture<Texture>("alpha", 0.1f);
        }

        if (props.has_property("specular_reflectance"))
            m_specular_reflectance = props.texture<Texture>("specular_reflectance", 1.f);

        m_flags = BSDFFlags::GlossyReflection | BSDFFlags::FrontSide;
        if (m_alpha_u != m_alpha_v)
            m_flags = m_flags | BSDFFlags::Anisotropic;

        m_components.clear();
        m_components.push_back(m_flags);
    }

    void traverse(TraversalCallback *callback) override {
        if (m_specular_reflectance)
            callback->put_object("specular_reflectance", m_specular_reflectance.get(), +ParamFlags::Differentiable);
        if (!has_flag(m_flags, BSDFFlags::Anisotropic))
            callback->put_object("alpha",                m_alpha_u.get(),              ParamFlags::Differentiable | ParamFlags::Discontinuous);
        else {
            callback->put_object("alpha_u",              m_alpha_u.get(),              ParamFlags::Differentiable | ParamFlags::Discontinuous);
            callback->put_object("alpha_v",              m_alpha_v.get(),              ParamFlags::Differentiable | ParamFlags::Discontinuous);
        }
        callback->put_object("eta", m_eta.get(), ParamFlags::Differentiable | ParamFlags::Discontinuous);
        callback->put_object("k",   m_k.get(),   ParamFlags::Differentiable | ParamFlags::Discontinuous);
    }

    // Jakob2014のプログラム-------------------------------------------------------------------------
    
    // ベクトルと行列の積を計算する関数
    dr::Array<float, 3> vecmatmul(dr::Array<float, 3> vec, dr::Matrix<float, 3> mat) {
        dr::Array<float, 3> result;
        for(int i = 0; i < 3; i++){
            for(int j = 0; j < 3; j++){
                result[i] = result[i] + vec[j] * mat[j][i];
            }
        }
        return result;
    }

    // 式(9)の行列Cを求める関数
    bool equation9(dr::Array<float, 3> wi, dr::Array<float, 3> wo,
                    float gamma, dr::Array<float, 3> m) {


        dr::Array<float, 3> x_hat = dr::normalize(dr::cross(wi, wo));
        dr::Array<float, 3> y_hat = dr::normalize(wi - wo);
        dr::Array<float, 3> z_hat = dr::normalize(wi + wo);

        // gammaは度数なのでラジアンに直す
        float g_cos = std::cosf(gamma * (M_PI / 180.0f));
        float lambda1 = (dr::dot(wi, wo) + g_cos) / (1.0f - g_cos);
        float lambda2 = powf(1 / std::tanf((gamma / 2) * (M_PI / 180.0f)), 2.0f);

        dr::Matrix<float, 3> Q = (x_hat[0], y_hat[0], z_hat[0],
                                x_hat[1], y_hat[1], z_hat[1],
                                x_hat[2], y_hat[2], z_hat[2]);
        

        dr::Matrix<float, 3> Lambda = (lambda1, 0, 0,
                                    0, lambda2, 0,
                                    0, 0, -1.0f);
        
        dr::Matrix<float, 3> C = Q * Lambda * dr::transpose(Q);
        
        float result = dr::dot(vecmatmul(m, C), m);

        if(result <= 0) {
            return true;
        } 
        else {
            return false;
        }
    }

    // ---------------------------------------------------------------------------------------------

    std::pair<BSDFSample3f, Spectrum> sample(const BSDFContext &ctx,
                                             const SurfaceInteraction3f &si,
                                             Float /* sample1 */,
                                             const Point2f &sample2,
                                             Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::BSDFSample, active);

        BSDFSample3f bs = dr::zeros<BSDFSample3f>();
        Float cos_theta_i = Frame3f::cos_theta(si.wi);
        active &= cos_theta_i > 0.f;

        if (unlikely(!ctx.is_enabled(BSDFFlags::GlossyReflection) || dr::none_or<false>(active)))
            return { bs, 0.f };

        /* Construct a microfacet distribution matching the
           roughness values at the current surface position. */
        MicrofacetDistribution distr(m_type,
                                     m_alpha_u->eval_1(si, active),
                                     m_alpha_v->eval_1(si, active),
                                     m_sample_visible);

        // Sample M, the microfacet normal
        Normal3f m;
        std::tie(m, bs.pdf) = distr.sample(si.wi, sample2);

        // Perfect specular reflection based on the microfacet normal
        bs.wo = reflect(si.wi, m);
        bs.eta = 1.f;
        bs.sampled_component = 0;
        bs.sampled_type = +BSDFFlags::GlossyReflection;

        // Ensure that this is a valid sample
        active &= dr::neq(bs.pdf, 0.f) && Frame3f::cos_theta(bs.wo) > 0.f;

        UnpolarizedSpectrum weight;
        if (likely(m_sample_visible))
            weight = distr.smith_g1(bs.wo, m);
        else
            weight = distr.G(si.wi, bs.wo, m) * dr::dot(si.wi, m) /
                     (cos_theta_i * Frame3f::cos_theta(m));

        // Jacobian of the half-direction mapping
        bs.pdf /= 4.f * dr::dot(bs.wo, m);

        // Evaluate the Fresnel factor
        dr::Complex<UnpolarizedSpectrum> eta_c(m_eta->eval(si, active),
                                           m_k->eval(si, active));

        Spectrum F;
        if constexpr (is_polarized_v<Spectrum>) {
            /* Due to the coordinate system rotations for polarization-aware
               pBSDFs below we need to know the propagation direction of light.
               In the following, light arrives along `-wo_hat` and leaves along
               `+wi_hat`. */
            Vector3f wo_hat = ctx.mode == TransportMode::Radiance ? bs.wo : si.wi,
                     wi_hat = ctx.mode == TransportMode::Radiance ? si.wi : bs.wo;

            // Mueller matrix for specular reflection.
            F = mueller::specular_reflection(UnpolarizedSpectrum(dot(wo_hat, m)), eta_c);

            /* The Stokes reference frame vector of this matrix lies perpendicular
               to the plane of reflection. */
            Vector3f s_axis_in  = dr::cross(m, -wo_hat);
            Vector3f s_axis_out = dr::cross(m, wi_hat);

            // Singularity when the input & output are collinear with the normal
            Mask collinear = dr::all(dr::eq(s_axis_in, Vector3f(0)));
            s_axis_in  = dr::select(collinear, Vector3f(1, 0, 0),
                                               dr::normalize(s_axis_in));
            s_axis_out = dr::select(collinear, Vector3f(1, 0, 0),
                                               dr::normalize(s_axis_out));

            /* Rotate in/out reference vector of F s.t. it aligns with the implicit
               Stokes bases of -wo_hat & wi_hat. */
            F = mueller::rotate_mueller_basis(F,
                                              -wo_hat, s_axis_in, mueller::stokes_basis(-wo_hat),
                                               wi_hat, s_axis_out, mueller::stokes_basis(wi_hat));
        } else {
            F = fresnel_conductor(UnpolarizedSpectrum(dr::dot(si.wi, m)), eta_c);
        }

        /* If requested, include the specular reflectance component */
        if (m_specular_reflectance)
            weight *= m_specular_reflectance->eval(si, active);

        return { bs, (F * weight) & active };
    }

    Spectrum eval(const BSDFContext &ctx, const SurfaceInteraction3f &si,
                  const Vector3f &wo, Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::BSDFEvaluate, active);

        Float cos_theta_i = Frame3f::cos_theta(si.wi),
              cos_theta_o = Frame3f::cos_theta(wo);

        active &= cos_theta_i > 0.f && cos_theta_o > 0.f;

        if (unlikely(!ctx.is_enabled(BSDFFlags::GlossyReflection) || dr::none_or<false>(active)))
            return 0.f;

        // Calculate the half-direction vector
        Vector3f H = dr::normalize(wo + si.wi);

        /* Construct a microfacet distribution matching the
           roughness values at the current surface position. */
        MicrofacetDistribution distr(m_type,
                                     m_alpha_u->eval_1(si, active),
                                     m_alpha_v->eval_1(si, active),
                                     m_sample_visible);

        // Evaluate the microfacet normal distribution
        Float D = distr.eval(H);

        active &= dr::neq(D, 0.f);

        // Evaluate Smith's shadow-masking function
        Float G = distr.G(si.wi, wo, H);

        // Evaluate the full microfacet model (except Fresnel)
        UnpolarizedSpectrum result = D * G / (4.f * Frame3f::cos_theta(si.wi));

        // Evaluate the Fresnel factor
        dr::Complex<UnpolarizedSpectrum> eta_c(m_eta->eval(si, active),
                                           m_k->eval(si, active));

        Spectrum F;
        if constexpr (is_polarized_v<Spectrum>) {
            /* Due to the coordinate system rotations for polarization-aware
               pBSDFs below we need to know the propagation direction of light.
               In the following, light arrives along `-wo_hat` and leaves along
               `+wi_hat`. */
            Vector3f wo_hat = ctx.mode == TransportMode::Radiance ? wo : si.wi,
                     wi_hat = ctx.mode == TransportMode::Radiance ? si.wi : wo;

            // Mueller matrix for specular reflection.
            F = mueller::specular_reflection(UnpolarizedSpectrum(dot(wo_hat, H)), eta_c);

            /* The Stokes reference frame vector of this matrix lies perpendicular
               to the plane of reflection. */
            Vector3f s_axis_in  = dr::cross(H, -wo_hat);
            Vector3f s_axis_out = dr::cross(H, wi_hat);

            // Singularity when the input & output are collinear with the normal
            Mask collinear = dr::all(dr::eq(s_axis_in, Vector3f(0)));
            s_axis_in  = dr::select(collinear, Vector3f(1, 0, 0),
                                               dr::normalize(s_axis_in));
            s_axis_out = dr::select(collinear, Vector3f(1, 0, 0),
                                               dr::normalize(s_axis_out));

            /* Rotate in/out reference vector of F s.t. it aligns with the implicit
               Stokes bases of -wo_hat & wi_hat. */
            F = mueller::rotate_mueller_basis(F,
                                              -wo_hat, s_axis_in, mueller::stokes_basis(-wo_hat),
                                               wi_hat, s_axis_out, mueller::stokes_basis(wi_hat));
        } else {
            F = fresnel_conductor(UnpolarizedSpectrum(dr::dot(si.wi, H)), eta_c);
        }

        /* If requested, include the specular reflectance component */
        if (m_specular_reflectance)
            result *= m_specular_reflectance->eval(si, active);

        return (F * result) & active;
    }

    Float pdf(const BSDFContext &ctx, const SurfaceInteraction3f &si,
              const Vector3f &wo, Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::BSDFEvaluate, active);

        Float cos_theta_i = Frame3f::cos_theta(si.wi),
              cos_theta_o = Frame3f::cos_theta(wo);

        // Calculate the half-direction vector
        Vector3f m = dr::normalize(wo + si.wi);

        /* Filter cases where the micro/macro-surface don't agree on the side.
           This logic is evaluated in smith_g1() called as part of the eval()
           and sample() methods and needs to be replicated in the probability
           density computation as well. */
        active &= cos_theta_i > 0.f && cos_theta_o > 0.f &&
                  dr::dot(si.wi, m) > 0.f && dr::dot(wo, m) > 0.f;

        if (unlikely(!ctx.is_enabled(BSDFFlags::GlossyReflection) || dr::none_or<false>(active)))
            return 0.f;

        /* Construct a microfacet distribution matching the
           roughness values at the current surface position. */
        MicrofacetDistribution distr(m_type,
                                     m_alpha_u->eval_1(si, active),
                                     m_alpha_v->eval_1(si, active),
                                     m_sample_visible);

        Float result;
        if (likely(m_sample_visible))
            result = distr.eval(m) * distr.smith_g1(si.wi, m) /
                     (4.f * cos_theta_i);
        else
            result = distr.pdf(si.wi, m) / (4.f * dr::dot(wo, m));

        return dr::select(active, result, 0.f);
    }

    std::pair<Spectrum, Float> eval_pdf(const BSDFContext &ctx,
                                        const SurfaceInteraction3f &si,
                                        const Vector3f &wo,
                                        Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::BSDFEvaluate, active);

        Float cos_theta_i = Frame3f::cos_theta(si.wi),
              cos_theta_o = Frame3f::cos_theta(wo);

        // Calculate the half-direction vector
        Vector3f H = dr::normalize(wo + si.wi);

        /* Filter cases where the micro/macro-surface don't agree on the side.
           This logic is evaluated in smith_g1() called as part of the eval()
           and sample() methods and needs to be replicated in the probability
           density computation as well. */
        active &= cos_theta_i > 0.f && cos_theta_o > 0.f &&
                  dr::dot(si.wi, H) > 0.f && dr::dot(wo, H) > 0.f;

        if (unlikely(!ctx.is_enabled(BSDFFlags::GlossyReflection) || dr::none_or<false>(active)))
            return { 0.f, 0.f };

        /* Construct a microfacet distribution matching the
           roughness values at the current surface position. */
        MicrofacetDistribution distr(m_type,
                                     m_alpha_u->eval_1(si, active),
                                     m_alpha_v->eval_1(si, active),
                                     m_sample_visible);

        // Evaluate the microfacet normal distribution
        Float D = distr.eval(H);

        active &= dr::neq(D, 0.f);

        // Evaluate Smith's shadow-masking function
        Float smith_g1_wi = distr.smith_g1(si.wi, H);
        Float G = smith_g1_wi * distr.smith_g1(wo, H);

        // Evaluate the full microfacet model (except Fresnel)
        UnpolarizedSpectrum value = D * G / (4.f * Frame3f::cos_theta(si.wi));

        // Evaluate the Fresnel factor
        dr::Complex<UnpolarizedSpectrum> eta_c(m_eta->eval(si, active),
                                           m_k->eval(si, active));

        Spectrum F;
        if constexpr (is_polarized_v<Spectrum>) {
            /* Due to the coordinate system rotations for polarization-aware
               pBSDFs below we need to know the propagation direction of light.
               In the following, light arrives along `-wo_hat` and leaves along
               `+wi_hat`. */
            Vector3f wo_hat = ctx.mode == TransportMode::Radiance ? wo : si.wi,
                     wi_hat = ctx.mode == TransportMode::Radiance ? si.wi : wo;

            // Mueller matrix for specular reflection.
            F = mueller::specular_reflection(UnpolarizedSpectrum(dot(wo_hat, H)), eta_c);

            /* The Stokes reference frame vector of this matrix lies perpendicular
               to the plane of reflection. */
            Vector3f s_axis_in  = dr::cross(H, -wo_hat);
            Vector3f s_axis_out = dr::cross(H, wi_hat);

            // Singularity when the input & output are collinear with the normal
            Mask collinear = dr::all(dr::eq(s_axis_in, Vector3f(0)));
            s_axis_in  = dr::select(collinear, Vector3f(1, 0, 0),
                                               dr::normalize(s_axis_in));
            s_axis_out = dr::select(collinear, Vector3f(1, 0, 0),
                                               dr::normalize(s_axis_out));

            /* Rotate in/out reference vector of F s.t. it aligns with the implicit
               Stokes bases of -wo_hat & wi_hat. */
            F = mueller::rotate_mueller_basis(F,
                                              -wo_hat, s_axis_in, mueller::stokes_basis(-wo_hat),
                                               wi_hat, s_axis_out, mueller::stokes_basis(wi_hat));
        } else {
            F = fresnel_conductor(UnpolarizedSpectrum(dr::dot(si.wi, H)), eta_c);
        }

        // If requested, include the specular reflectance component
        if (m_specular_reflectance)
            value *= m_specular_reflectance->eval(si, active);

        Float pdf;
        if (likely(m_sample_visible))
            pdf = D * smith_g1_wi / (4.f * cos_theta_i);
        else
            pdf = distr.pdf(si.wi, H) / (4.f * dr::dot(wo, H));

        return { F * value & active, dr::select(active, pdf, 0.f) };
    }

    std::string to_string() const override {
        std::ostringstream oss;
        oss << "DiscreteRoughConductor[" << std::endl
            << "  distribution = " << m_type << "," << std::endl
            << "  sample_visible = " << m_sample_visible << "," << std::endl
            << "  alpha_u = " << string::indent(m_alpha_u) << "," << std::endl
            << "  alpha_v = " << string::indent(m_alpha_v) << "," << std::endl;
        if (m_specular_reflectance)
           oss << "  specular_reflectance = " << string::indent(m_specular_reflectance) << "," << std::endl;
        oss << "  eta = " << string::indent(m_eta) << "," << std::endl
            << "  k = " << string::indent(m_k) << std::endl
            << "]";
        return oss.str();
    }

    MI_DECLARE_CLASS()
private:
    /// Specifies the type of microfacet distribution
    MicrofacetType m_type;
    /// Anisotropic roughness values
    ref<Texture> m_alpha_u, m_alpha_v;
    /// Importance sample the distribution of visible normals?
    bool m_sample_visible;
    /// Relative refractive index (real component)
    ref<Texture> m_eta;
    /// Relative refractive index (imaginary component).
    ref<Texture> m_k;
    /// Specular reflectance component
    ref<Texture> m_specular_reflectance;
};

MI_IMPLEMENT_CLASS_VARIANT(DiscreteRoughConductor, BSDF)
MI_EXPORT_PLUGIN(DiscreteRoughConductor, "Discrete rough conductor")
NAMESPACE_END(mitsuba)
