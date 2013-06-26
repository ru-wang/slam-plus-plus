/*
								+--------------------------------+
								|                                |
								|       ***   Bitmap   ***       |
								|                                |
								| Copyright (c) -tHE SWINe- 2011 |
								|                                |
								|            Bitmap.h            |
								|                                |
								+--------------------------------+
*/

#pragma once
#ifndef __BITMAP_STRUCTURE_INCLUDED
#define __BITMAP_STRUCTURE_INCLUDED

/**
 *	@file Bitmap.h
 *	@date 2011
 *	@author -tHE SWINe-
 *	@brief a simple, easy to use bitmap class
 *
 *	@todo allow for use of shaders (textures) for line and triangle rasterization (even without z-buffer)
 *	@todo write code for rendering axis-aligned rectangles (circles?)
 *	@todo write code for rendering antialiassed lines
 *
 *	@date 2012-06-19
 *
 *	Added \#pragma once.
 *
 *	@date 2012-11-02
 *
 *	Added subpixel precise and antialiassed line rasterization routines.
 *	Added axis aligned rectangle rasterization routines.
 *
 *	@note All the new functionality is largerly untested - test it.
 *
 */

/**
 *	@def __BMP_INCLUDED
 *	@brief legacy header guard name
 */
#define __BMP_INCLUDED

#include <new> // nothrow_t
#include <utility> // swap
#include <algorithm> // std::min, std::max
#include <string.h> // memcpy
#include <math.h>
#include "slam/Integer.h"

/**
 *	@brief a simple raster image with fixed RGBA8 storage
 *
 *	The image data are always RGBA8. Alpha is stored in the most significant byte,
 *	followed by red, green and blue with decreasing significance.
 *
 *	The storage is very simple, each 32 bits in the buffer contains a single pixel,
 *	the first pixel is in top left corner, there is no scanline padding.
 */
struct TBmp {
	char n_former_bpp; /**< @brief former bpp, before conversion to RGBA8 */
	bool b_grayscale; /**< @brief grayscale flag (if set, the bitmap is assumed
		to contain grayscale image, stored as RGBA8) */
	bool b_alpha; /**< @brief alpha channel flag (if set, the alpha channel is significant;
		otherwise it's expected to be 0xff in all image pixels) */
	int n_width; /**< @brief image width, in pixels */
	int n_height; /**< @brief image height, in pixels */
	uint32_t *p_buffer; /**< @brief pointer to image data */

	/**
	 *	@brief deletes data buffer
	 */
	inline void Free()
	{
		delete[] p_buffer;
	}

	/**
	 *	@brief deletes data buffer and this
	 */
	inline void Delete()
	{
		delete[] p_buffer;
		delete this;
	}

	/**
	 *	@brief allocates a new bitmap with the specified parameters
	 *
	 *	@param[in] _n_width is width of the bitmap, in pixels
	 *	@param[in] _n_height is height of the bitmap, in pixels
	 *	@param[in] _b_grayscale is grayscale flag
	 *	@param[in] _b_alpha is alpha channel flag
	 *	@param[in] _n_former_bpp is number of bits of the intended contents
	 *		of the bitmap, before conversion to RGBA8
	 *
	 *	@return Returns pointer to the new bitmap on success, or 0 on failure.
	 */
	static TBmp *p_Alloc(int _n_width, int _n_height,
		bool _b_grayscale = false, bool _b_alpha = false, int _n_former_bpp = 8)
	{
		TBmp *p_bmp;
		if(!(p_bmp = new(std::nothrow) TBmp))
			return 0;
		p_bmp->n_width = _n_width;
		p_bmp->n_height = _n_height;
		p_bmp->b_grayscale = _b_grayscale;
		p_bmp->b_alpha = _b_alpha;
		p_bmp->n_former_bpp = _n_former_bpp;
		if(!(p_bmp->p_buffer = new(std::nothrow) uint32_t[size_t(_n_width) * _n_height])) {
			delete p_bmp;
			return 0;
		}
		return p_bmp;
	}

	/**
	 *	@copydoc p_Alloc
	 *	@deprecated This is deprecated in favor of p_Alloc(), the effect is the same.
	 */
	static inline TBmp *p_CreateBitmap(int _n_width, int _n_height,
		bool _b_grayscale = false, bool _b_alpha = false, int _n_former_bpp = 8)
	{
		return p_Alloc(_n_width, _n_height, _b_grayscale, _b_alpha, _n_former_bpp);
	}

	/**
	 *	@brief creates a clone of the image
	 *	@param[in] b_ignore_contents is ignore contents flag (if set, only bitmap with
	 *		the same parameters is created, but the raster data are not copied)
	 *	@return Returns pointer to the new bitmap on success, or 0 on failure.
	 */
	TBmp *p_Clone(bool b_ignore_contents = false) const
	{
		TBmp *p_clone = p_CreateBitmap(n_width, n_height, b_grayscale, b_alpha, n_former_bpp);
		if(p_buffer && p_clone && !b_ignore_contents)
			memcpy(p_clone->p_buffer, p_buffer, size_t(n_width) * n_height * sizeof(uint32_t));
		return p_clone;
	}

	/**
	 *	@brief crops the image
	 *
	 *	@param[in] x is a coordinate of top-left corner of the crop rectangle
	 *	@param[in] y is a coordinate of top-left corner of the crop rectangle
	 *	@param[in] _n_width is width of the crop rectangle, in pixels
	 *	@param[in] _n_height is height of the crop rectangle, in pixels
	 *
	 *	@return Returns pointer to the new bitmap on success, or 0 on failure.
	 *
	 *	@note The cropping rectangle must lie completely inside the bitmap.
	 */
	TBmp *p_Crop(int x, int y, int _n_width, int _n_height) const
	{
		_ASSERTE(x >= 0 && y >= 0 && x + _n_width <= n_width && y + _n_height <= n_height);

		TBmp *p_clone = p_CreateBitmap(_n_width, _n_height, b_grayscale, b_alpha, n_former_bpp);
		if(p_buffer && p_clone) {
			for(int yy = 0; yy < _n_height; ++ yy) {
				memcpy(p_clone->p_buffer + _n_width * yy,
					p_buffer + (x + (y + yy) * n_width), _n_width * sizeof(uint32_t));
			}
		}
		return p_clone;
	}

	/**
	 *	@brief fills the bitmap with constant color
	 *	@param[in] n_color is the fill color
	 */
	void Clear(uint32_t n_color)
	{
		for(uint32_t *p_dest = p_buffer, *p_end = p_buffer + size_t(n_width) * n_height; p_dest != p_end; ++ p_dest)
			*p_dest = n_color;
	}

	/**
	 *	@brief converts the bitmap to grayscale
	 *	@note This has no effect if the b_grayscale flag is already set.
	 */
	void Make_Grayscale()
	{
		if(b_grayscale)
			return;
		// that was easy ...

		for(uint32_t *p_pixel = p_buffer, *p_end = p_buffer + (size_t(n_width) * n_height); p_pixel != p_end; ++ p_pixel) {
			uint32_t c = *p_pixel;
			int r = (c >> 16) & 0xff;
			int g = (c >> 8) & 0xff;
			int b = c & 0xff;
			int n_grey = (r * int(.299f * 0x10000) + g * int(.578f * 0x10000) + b * int(.114f * 0x10000)) >> 16;
			*p_pixel = 0xff000000U | (n_grey << 16) | (n_grey << 8) | n_grey;
		}
		b_grayscale = true;
	}

	/**
	 *	@brief fills a selected pixel with a given color
	 *
	 *	@param[in] x is a coordinate of the pixel to fill
	 *	@param[in] y is a coordinate of the pixel to fill
	 *	@param[in] n_color is the fill color
	 *
	 *	@note This performs array boundary checking,
	 *		coordinates outside the bitmap are ok.
	 */
	inline void PutPixel(int x, int y, uint32_t n_color)
	{
		if(x >= 0 && x < n_width && y >= 0 && y < n_height)
			p_buffer[int(x) + n_width * int(y)] = n_color;
	}

	/**
	 *	@brief draws an axis aligned rectangle (only lines, no fill)
	 *
	 *	@param[in] n_x0 is a coordinate of the top-left corner
	 *	@param[in] n_y0 is a coordinate of the top-left corner
	 *	@param[in] n_x1 is a coordinate of the bottom-right corner
	 *	@param[in] n_y1 is a coordinate of the bottom-right corner
	 *	@param[in] n_color is the line color
	 *
	 *	@todo Add support for line width.
	 */
	void DrawRect(int n_x0, int n_y0, int n_x1, int n_y1, uint32_t n_color)
	{
		if(n_x0 > n_x1)
			std::swap(n_x0, n_x1);
		if(n_y0 > n_y1)
			std::swap(n_y0, n_y1);
		// make sure it is ordered

		if(n_y0 >= 0 && n_y0 < n_height) {
			_ASSERTE(n_y1 >= 0);
			if(n_y1 < n_height) {
				for(int x = std::max(0, n_x0); x < std::min(n_width - 1, n_x1); ++ x) {
					p_buffer[x + n_y0 * n_width] = n_color;
					p_buffer[x + n_y1 * n_width] = n_color;
				}
				// both are in
			} else {
				for(int x = std::max(0, n_x0); x < std::min(n_width - 1, n_x1); ++ x)
					p_buffer[x + n_y0 * n_width] = n_color;
			}
		} else if(n_y1 >= 0 && n_y1 < n_height) {
			for(int x = std::max(0, n_x0); x < std::min(n_width - 1, n_x1); ++ x)
				p_buffer[x + n_y1 * n_width] = n_color;
		}
		// draw horizontal lines

		if(n_x0 >= 0 && n_x0 < n_width) {
			_ASSERTE(n_x1 >= 0);
			if(n_x1 < n_width) {
				for(int y = std::max(0, n_y0); y < std::min(n_height - 1, n_y1); ++ y) {
					p_buffer[n_x0 + y * n_width] = n_color;
					p_buffer[n_x1 + y * n_width] = n_color;
				}
				// both are in
			} else {
				for(int y = std::max(0, n_y0); y < std::min(n_height - 1, n_y1); ++ y)
					p_buffer[n_x0 + y * n_width] = n_color;
			}
		} else if(n_x1 >= 0 && n_x1 < n_width) {
			for(int y = std::max(0, n_y0); y < std::min(n_height - 1, n_y1); ++ y)
				p_buffer[n_x1 + y * n_width] = n_color;
		}
		// draw vertical lines
	}

	/**
	 *	@brief fills an axis aligned rectangle with constant color
	 *
	 *	@param[in] n_x0 is a coordinate of the top-left corner
	 *	@param[in] n_y0 is a coordinate of the top-left corner
	 *	@param[in] n_x1 is a coordinate of the bottom-right corner
	 *	@param[in] n_y1 is a coordinate of the bottom-right corner
	 *	@param[in] n_color is the fill color
	 */
	void FillRect(int n_x0, int n_y0, int n_x1, int n_y1, uint32_t n_color)
	{
		if(n_x0 > n_x1)
			std::swap(n_x0, n_x1);
		if(n_y0 > n_y1)
			std::swap(n_y0, n_y1);
		// make sure it is ordered

		if(n_x1 < 0 || n_y1 < 0)
			return;
		if(n_x0 >= n_width || n_y0 >= n_height)
			return;
		// simple rejection

		n_x0 = std::max(0, n_x0);
		n_y0 = std::max(0, n_y0);
		n_x1 = std::min(n_width, n_x1 + 1) - n_x0; // number of pixels to fill
		_ASSERTE(n_x1 >= 0);
		_ASSERTE(n_x0 + n_x1 <= n_width);
		n_y1 = std::min(n_height, n_y1 + 1); // y one past the last scanline
		// make sure it is inside

		uint32_t *p_scan = p_buffer + n_x0 + n_y0 * n_width;
		for(int y = n_y0; y < n_y1; ++ y, p_scan += n_width) {
			for(uint32_t *p_ptr = p_scan, *p_end = p_scan + n_x1; p_ptr != p_end; ++ p_ptr)
				*p_ptr = n_color;
		}
		// fill
	}

	/**
	 *	@brief clips a line to bitmap interior
	 *
	 *	@param[in,out] r_f_x0 is a coordinate of the first line point
	 *	@param[in,out] r_f_y0 is a coordinate of the first line point
	 *	@param[in,out] r_f_x1 is a coordinate of the second line point
	 *	@param[in,out] r_f_y1 is a coordinate of the second line point
	 *
	 *	@return Returns true if the line is inside, false if it is completely
	 *		outside (early reject, values of arguments are not changed).
	 */
	inline bool ClipLine(float &r_f_x0, float &r_f_y0, float &r_f_x1, float &r_f_y1) const
	{
		bool b_not_narrow;
		float f_dxdy = ((b_not_narrow = (fabs(r_f_x1 - r_f_x0) > 1e-5f)))?
			(r_f_y1 - r_f_y0) / (r_f_x1 - r_f_x0) : 0;
		if(r_f_x0 < 0 || r_f_x1 < 0) {
			if(r_f_x0 < 0 && r_f_x1 < 0)
				return false; // offscreen
			if(r_f_x0 < 0) {
				r_f_y0 -= f_dxdy * r_f_x0; // note this rounds ...
				r_f_x0 = 0;
			} else {
				r_f_y1 -= f_dxdy * r_f_x1; // note this rounds ...
				r_f_x1 = 0;
			}
		}
		const int n_w_max = n_width - 1;
		if(r_f_x0 > n_w_max || r_f_x1 > n_w_max) {
			if(r_f_x0 > n_w_max && r_f_x1 > n_w_max)
				return false; // offscreen
			if(r_f_x0 > n_w_max) {
				float dx = r_f_x0 - n_w_max;
				r_f_y0 -= f_dxdy * dx; // note this rounds ...
				r_f_x0 = float(n_w_max);
			} else {
				float dx = r_f_x1 - n_w_max;
				r_f_y1 -= f_dxdy * dx; // note this rounds ...
				r_f_x1 = float(n_w_max);
			}
		}
		if(!b_not_narrow)
			f_dxdy = 1e37f; // stable value for this part (or could branch below)
		if(r_f_y0 < 0 || r_f_y1 < 0) {
			if(r_f_y0 < 0 && r_f_y1 < 0)
				return false; // offscreen
			if(r_f_y0 < 0) {
				r_f_x0 -= r_f_y0 / f_dxdy; // note this rounds ...
				r_f_y0 = 0;
			} else {
				r_f_x1 -= r_f_y1 / f_dxdy; // note this rounds ...
				r_f_y1 = 0;
			}
		}
		const int n_h_max = n_height - 1;
		if(r_f_y0 > n_h_max || r_f_y1 > n_h_max) {
			if(r_f_y0 > n_h_max && r_f_y1 > n_h_max)
				return false; // offscreen
			if(r_f_y0 > n_h_max) {
				float dy = r_f_y0 - n_h_max;
				r_f_x0 -= dy / f_dxdy; // note this rounds ...
				r_f_y0 = float(n_h_max);
			} else {
				float dy = r_f_y1 - n_h_max;
				r_f_x1 -= dy / f_dxdy; // note this rounds ...
				r_f_y1 = float(n_h_max);
			}
		}
		// perform simple clipping

		_ASSERTE(int(r_f_x0) >= 0 && int(r_f_x0) <= n_w_max);
		_ASSERTE(int(r_f_y0) >= 0 && int(r_f_y0) <= n_h_max);
		_ASSERTE(int(r_f_x1) >= 0 && int(r_f_x1) <= n_w_max);
		_ASSERTE(int(r_f_y1) >= 0 && int(r_f_y1) <= n_h_max);

		return true;
	}

	/**
	 *	@brief a simple interpolator object for rasterizing lines
	 */
	class CInterpolator {
	public:
		typedef int TFixedPoint; /**< @brief this interpolator should use fixed-point numbers for increased accuracy */

	protected:
		TFixedPoint m_y; /**< @brief current interpolated value */
		int m_n_length; /**< @brief length of the interpolation domain */
		TFixedPoint m_n_slope; /**< @brief interpolated line slope */
		TFixedPoint m_n_error; /**< @brief error accumulator */
		TFixedPoint m_n_remainder; /**< @brief error increase per interpolation step */

	public:
		/**
		 *	@brief default constructor
		 *
		 *	@param[in] n_y0 is the interpolated value at the beginning of the domain
		 *	@param[in] n_y1 is the interpolated value at the end of the domain
		 *	@param[in] n_length is length of the domain
		 */
		inline CInterpolator(TFixedPoint n_y0, TFixedPoint n_y1, int n_length)
			:m_y(n_y0), m_n_length(std::max(1, n_length)), m_n_slope((n_y1 - n_y0) / m_n_length),
			m_n_error((n_y1 - n_y0) % m_n_length), m_n_remainder(m_n_error)
		{
			if(m_n_error <= 0) {
				m_n_error += n_length;
				m_n_remainder += n_length;
				-- m_n_slope;
			}
			// fix the rounding direction for negative values

			m_n_error -= n_length; // !!
		}

		/**
		 *	@brief calculates interpolated value at the next position
		 */
		inline void operator ++()
		{
			m_y += m_n_slope;
			if((m_n_error += m_n_remainder) > 0) {
				m_n_error -= m_n_length;
				++ m_y;
			}
			// DDA
		}

		/**
		 *	@brief gets interpolated value
		 *	@return Returns the current interpolated value.
		 */
		inline TFixedPoint y() const
		{
			return m_y;
		}
	};

	/**
	 *	@brief draws a solid-color line
	 *
	 *	@param[in] n_x0 is a coordinate of the first line point
	 *	@param[in] n_y0 is a coordinate of the first line point
	 *	@param[in] n_x1 is a coordinate of the second line point
	 *	@param[in] n_y1 is a coordinate of the second line point
	 *	@param[in] n_color is the line color
	 *	@param[in] n_line_width is the line width, in pixels
	 *
	 *	@note This performs clipping, coordinates outside the bitmap are ok.
	 */
	void DrawLine(int n_x0, int n_y0, int n_x1, int n_y1, uint32_t n_color, int n_line_width = 1)
	{
		if(n_line_width <= 0)
			return;
		// too thin

		float f_x0 = float(n_x0) + .5f, f_y0 = float(n_y0) + .5f;
		float f_x1 = float(n_x1) + .5f, f_y1 = float(n_y1) + .5f;
		// adjust lines with integer coordinates to be on pixel centers (makes tilted lines look better)

		DrawLine_SP(f_x0, f_y0, f_x1, f_y1, n_color, n_line_width);
		// use fancy function
	}

	/**
	 *	@brief draws a solid-color line with subpixel precision
	 *
	 *	@param[in] f_x0 is a coordinate of the first line point
	 *	@param[in] f_y0 is a coordinate of the first line point
	 *	@param[in] f_x1 is a coordinate of the second line point
	 *	@param[in] f_y1 is a coordinate of the second line point
	 *	@param[in] n_color is the line color
	 *	@param[in] n_line_width is the line width, in pixels
	 *
	 *	@note This performs clipping, coordinates outside the bitmap are ok.
	 */
	void DrawLine_SP(float f_x0, float f_y0, float f_x1, float f_y1,
		uint32_t n_color, int n_line_width = 1)
	{
		if(n_line_width <= 0)
			return;
		// too thin

		if(!ClipLine(f_x0, f_y0, f_x1, f_y1))
			return;
		// perform simple clipping

		/*if(f_x0 >= 0 && f_x0 < n_width && f_y0 >= 0 && f_y0 < n_height)
			p_buffer[int(f_x0) + n_width * int(f_y0)] = 0xffff00ff;
		if(f_x1 >= 0 && f_x1 < n_width && f_y1 >= 0 && f_y1 < n_height)
			p_buffer[int(f_x1) + n_width * int(f_y1)] = 0xffff00ff;*/
		// debug - mark endpoints // now it fails olny on very short lines (the test with the spiral)

		_ASSERTE(std::max(abs(int(f_x0) - int(f_x1)), abs(int(f_y0) - int(f_y1))) <=
			std::max(n_width, n_height)); // line lenght is now bound by bitmap size
		bool b_steep = fabs(f_y0 - f_y1) > fabs(f_x0 - f_x1);
		if(b_steep) {
			std::swap(f_x0, f_y0);
			std::swap(f_x1, f_y1);
			// makes sure it is rasterized in the larger dimension
		}
		if(f_x0 > f_x1) {
			std::swap(f_x0, f_x1);
			std::swap(f_y0, f_y1);
		}
		if(n_line_width == 1) {
			int n_len = abs(int(floor(f_x1) - floor(f_x0)));
			CInterpolator lerp(int(f_y0 * 256), int(f_y1 * 256), n_len);
			// note this is not adjusted for fractional x

			_ASSERTE(f_x0 <= f_x1);
			if(b_steep) {
				for(int n_y = int(floor(f_x0)), n_end = int(floor(f_x1)) + 1; n_y < n_end; ++ n_y, ++ lerp) {
					int n_x = lerp.y() >> 8;
					_ASSERTE(n_x >= 0 && n_x < n_width && n_y >= 0 && n_y < n_height);
					p_buffer[n_x + n_y * n_width] = n_color;
				}
			} else {
				for(int n_x = int(floor(f_x0)), n_end = int(floor(f_x1)) + 1; n_x < n_end; ++ n_x, ++ lerp) {
					int n_y = lerp.y() >> 8;
					_ASSERTE(n_x >= 0 && n_x < n_width && n_y >= 0 && n_y < n_height);
					p_buffer[n_x + n_y * n_width] = n_color;
				}
			}
			// thin lines
		} else {
			float f_dx = fabs(f_x1 - f_x0);
			float f_dy = fabs(f_y1 - f_y0);
			float f_thickness_scale = sqrt(f_dx * f_dx + f_dy * f_dy) / std::max(1.0f, std::max(f_dx, f_dy));
			n_line_width = int(n_line_width * f_thickness_scale + .5f);
			// adjust line thickness based on line angle

			int n_line_extent_top = (n_line_width - 1) / 2;
			int n_line_extent_bottom = (n_line_width - 1) - n_line_extent_top;
			// calculate extent on top and bottom

			int n_len = abs(int(floor(f_x1) - floor(f_x0)));
			CInterpolator lerp(int(f_y0 * 256), int(f_y1 * 256), n_len);
			// note this is not adjusted for fractional n_x

			if(b_steep) {
				for(int n_y = int(floor(f_x0)), n_end = int(floor(f_x1)) + 1; n_y < n_end; ++ n_y, ++ lerp) {
					int n_x = lerp.y() >> 8;
					_ASSERTE(n_x >= 0 && n_x < n_width && n_y >= 0 && n_y < n_height);
					p_buffer[n_x + n_y * n_width] = n_color;
					for(int dy = 1; dy <= n_line_extent_top; ++ dy) {
						if(n_x + dy < n_width)
							p_buffer[n_x + dy + n_y * n_width] = n_color;
					}
					for(int dy = 1; dy <= n_line_extent_bottom; ++ dy) {
						if(n_x >= dy)
							p_buffer[n_x - dy + n_y * n_width] = n_color;
					}
				}
			} else {
				for(int n_x = int(floor(f_x0)), n_end = int(floor(f_x1)) + 1; n_x < n_end; ++ n_x, ++ lerp) {
					int n_y = lerp.y() >> 8;
					_ASSERTE(n_x >= 0 && n_x < n_width && n_y >= 0 && n_y < n_height);
					p_buffer[n_x + n_y * n_width] = n_color;
					for(int dy = 1; dy <= n_line_extent_top; ++ dy) {
						if(n_y + dy < n_height)
							p_buffer[n_x + (n_y + dy) * n_width] = n_color;
					}
					for(int dy = 1; dy <= n_line_extent_bottom; ++ dy) {
						if(n_y >= dy)
							p_buffer[n_x + (n_y - dy) * n_width] = n_color;
					}
				}
			}
			// thick lines
		}

		/*if(b_steep) { // !!
			std::swap(f_x0, f_y0);
			std::swap(f_x1, f_y1);
		}
		if(f_x0 >= 0 && f_x0 < n_width && f_y0 >= 0 && f_y0 < n_height)
			if(p_buffer[int(f_x0) + n_width * int(f_y0)] != n_color) printf("start\n");
		if(f_x1 >= 0 && f_x1 < n_width && f_y1 >= 0 && f_y1 < n_height)
			if(p_buffer[int(f_x1) + n_width * int(f_y1)] != n_color) printf("end\n");*/
		// make sure that the endpoints were filled
	}

	/**
	 *	@brief modulates RGBA color by alpha
	 *
	 *	@param[in] n_src is a RGBA color
	 *	@param[in] n_alpha is modulation coefficient (in 0 to 255 range)
	 *
	 *	@return Returns RGBA color, modulated by alpha.
	 */
	static inline uint32_t n_Modulate(uint32_t n_src, int n_alpha)
	{
		_ASSERTE(n_alpha >= 0 && n_alpha <= 0xff);
		return (((((n_src & 0xff00ff) * n_alpha) & 0xff00ff00U) >> 8) |
			   ((((n_src & 0xff00ff00) >> 8) * n_alpha) & 0xff00ff00U));
		// use two ops on pairs of elems
	}

	/**
	 *	@brief modulates RGB color by alpha
	 *
	 *	@param[in] n_src is a RGB color (alpha is ignored)
	 *	@param[in] n_alpha is modulation coefficient (in 0 to 255 range)
	 *
	 *	@return Returns RGB color (alpha is null), modulated by alpha.
	 */
	static inline uint32_t n_Modulate_RGB(uint32_t n_src, int n_alpha)
	{
		_ASSERTE(n_alpha >= 0 && n_alpha <= 0xff);
		return ((((n_src & 0xff00ff) * n_alpha) & 0xff00ff00U) |
			   (((n_src & 0x00ff00) * n_alpha) & 0x00ff0000U)) >> 8; // alpha overflows, and is dammaged
		// use two ops on pairs of elems, save one bit shift
	}

	/**
	 *	@brief modulates R color by alpha
	 *
	 *	@param[in] n_src is a R color (GBA is ignored)
	 *	@param[in] n_alpha is modulation coefficient (in 0 to 255 range)
	 *
	 *	@return Returns R color (the other components are null), modulated by alpha.
	 */
	static inline uint32_t n_Modulate_Red(uint32_t n_src, int n_alpha)
	{
		_ASSERTE(n_alpha >= 0 && n_alpha <= 0xff);
		return (((n_src & 0xff) * n_alpha) & 0xff00) >> 8; // only the red channel is returned
	}

	/**
	 *	@brief modulates greyscale color by alpha
	 *
	 *	@param[in] n_src is a R color (GBA is ignored)
	 *	@param[in] n_alpha is modulation coefficient (in 0 to 255 range)
	 *
	 *	@return Returns greyscale color (RGB is red * alpha, A is alpha), modulated by alpha.
	 */
	static inline uint32_t n_Modulate_Grey(uint8_t n_src, int n_alpha)
	{
		_ASSERTE(n_alpha >= 0 && n_alpha <= 0xff);
		uint32_t n_grey = (uint32_t(n_src) * n_alpha) & 0xff00;
		return (n_alpha << 24) | n_grey | (n_grey << 8) | (n_grey >> 8); // only the red channel is returned
	}

	/**
	 *	@brief blends two RGBA colors based on alpha
	 *
	 *	Calculates r_n_dest = r_n_dest * (255 - n_alpha) + n_src * n_alpha, with RGBA arithmetic.
	 *
	 *	@param[in,out] r_n_dest is destination and left operand (the framebuffer)
	 *	@param[in] n_src is RGBA color (right operand)
	 *	@param[in] n_alpha is modulation coefficient (in 0 to 255 range)
	 */
	static inline void AlphaBlend(uint32_t &r_n_dest, uint32_t n_src, int n_alpha)
	{
		_ASSERTE(n_alpha >= 0 && n_alpha <= 0xff);
		r_n_dest = n_Modulate(r_n_dest, 0xff - n_alpha) + n_Modulate(n_src, n_alpha);
	}

	/**
	 *	@brief draws a solid-color antialiassed line with subpixel precision
	 *
	 *	@param[in] f_x0 is a coordinate of the first line point
	 *	@param[in] f_y0 is a coordinate of the first line point
	 *	@param[in] f_x1 is a coordinate of the second line point
	 *	@param[in] f_y1 is a coordinate of the second line point
	 *	@param[in] n_color is the line color
	 *	@param[in] n_line_width is the line width, in pixels
	 *
	 *	@note This performs clipping, coordinates outside the bitmap are ok.
	 *	@note This is somewhat wasteful if the bitmap is grayscale
	 *		since the blending full is RGBA blending.
	 *
	 *	@todo For thick line rasterization, the line endpoints are not finished,
	 *		seams might occur where lines with different slopes meet.
	 *	@todo Thick line rasterization is limited to integer line widths,
	 *		some widths (e.g. 2) do not give nice results.
	 */
	void DrawLine_AA(float f_x0, float f_y0, float f_x1, float f_y1,
		uint32_t n_color, int n_line_width = 1)
	{
		if(n_line_width <= 0)
			return;
		// too thin

		if(!ClipLine(f_x0, f_y0, f_x1, f_y1))
			return;
		// perform simple clipping

		_ASSERTE(std::max(abs(int(f_x0) - int(f_x1)), abs(int(f_y0) - int(f_y1))) <=
			std::max(n_width, n_height)); // line lenght is now bound by bitmap size
		bool b_steep = fabs(f_y0 - f_y1) > fabs(f_x0 - f_x1);
		if(b_steep) {
			std::swap(f_x0, f_y0);
			std::swap(f_x1, f_y1);
			// makes sure it is rasterized in the larger dimension
		}
		if(f_x0 > f_x1) {
			std::swap(f_x0, f_x1);
			std::swap(f_y0, f_y1);
		}
		float f_dxdy = (fabs(f_x1 - f_x0) > 1e-5f)? (f_y1 - f_y0) / (f_x1 - f_x0) : 0;
		int n_gradient = int(256 * f_dxdy);
		int n_end_x0 = int(floor(f_x0 + .5f));
		int n_end_x1 = int(floor(f_x1 + .5f));
		// note the .5 are important otherwise antialiassing discontinuities occur in the first quadrant

		if(n_line_width == 1) {
			if(n_end_x0 == n_end_x1) {
				float f_coverage = f_x1 - f_x0; // length of the line inside pixel
				float f_y_end = f_y0 + f_dxdy * (n_end_x0 - (f_x0 + f_x1) * .5f); // y-position of line center
				// average y in pixel

				int n_y_alpha = int(255 * (f_y_end - floor(f_y_end)));

				if(b_steep) {
					int n_y = n_end_x0, n_x = int(floor(f_y_end));
					if(n_x >= 0 && n_x < n_width && n_y >= 0 && n_y < n_height)
						AlphaBlend(p_buffer[n_x + n_y * n_width], n_color, int((255 - n_y_alpha) * f_coverage));
					if(n_x + 1 >= 0 && n_x + 1 < n_width && n_y >= 0 && n_y < n_height)
						AlphaBlend(p_buffer[n_x + 1 + n_y * n_width], n_color, int(n_y_alpha * f_coverage));
				} else {
					int n_x = n_end_x0, n_y = int(floor(f_y_end));
					if(n_x >= 0 && n_x < n_width && n_y >= 0 && n_y < n_height)
						AlphaBlend(p_buffer[n_x + n_y * n_width], n_color, int((255 - n_y_alpha) * f_coverage));
					if(n_x >= 0 && n_x < n_width && n_y + 1 >= 0 && n_y + 1 < n_height)
						AlphaBlend(p_buffer[n_x + (n_y + 1) * n_width], n_color, int(n_y_alpha * f_coverage));
				}

				return;
			}
			// in case the line only occupies a single pixel

			float f_end_y0 = f_y0 + f_dxdy * (n_end_x0 - f_x0);
			float f_cov_x0 = ceil(f_x0 + .5f) - (f_x0 + .5f); // how much of line is in the first pixel
			int n_lerp_y = int((f_end_y0 + f_dxdy) * 256);
			float f_end_y1 = f_y1 + f_dxdy * (n_end_x1 - f_x1);
			float f_cov_x1 = 1 - (ceil(f_x1 + .5f) - (f_x1 + .5f)); // how much of line is in the last pixel
			int n_alpha_y0 = 255 - int(255 * (ceil(f_end_y0) - f_end_y0));
			int n_alpha_y1 = int(255 * (f_end_y1 - floor(f_end_y1)));
			// calculate aliassing on the end of the lines
			// note the .5 are important otherwise antialiassing discontinuities occur in the first quadrant

			if(b_steep) {
				if(n_end_x0 >= 0 && n_end_x0 < n_height) {
					int n_y = n_end_x0, n_x = int(floor(f_end_y0));
					if(n_x >= 0 && n_x < n_width)
						AlphaBlend(p_buffer[n_x + n_y * n_width], n_color, int((255 - n_alpha_y0) * f_cov_x0));
					if(n_x + 1 >= 0 && n_x + 1 < n_width)
						AlphaBlend(p_buffer[n_x + 1 + n_y * n_width], n_color, int(n_alpha_y0 * f_cov_x0));
				}
				// handle the first endpoint

				if(n_end_x1 >= 0 && n_end_x1 < n_height) {
					int n_y = n_end_x1, n_x = int(floor(f_end_y1));
					if(n_x >= 0 && n_x < n_width)
						AlphaBlend(p_buffer[n_x + n_y * n_width], n_color, int((255 - n_alpha_y1) * f_cov_x1));
					if(n_x + 1 >= 0 && n_x + 1 < n_width)
						AlphaBlend(p_buffer[n_x + 1 + n_y * n_width], n_color, int(n_alpha_y1 * f_cov_x1));
				}
				// handle the second endpoint

				for(int n_y = n_end_x0 + 1; n_y < n_end_x1; ++ n_y) {
					int n_x = n_lerp_y >> 8;
					_ASSERTE(n_x >= 0 && n_x < n_width && n_y >= 0 && n_y < n_height);
					int n_alpha_0 = n_lerp_y & 0xff;
					AlphaBlend(p_buffer[n_x + n_y * n_width], n_color, 0xff - n_alpha_0);
					if(n_x + 1 < n_width)
						AlphaBlend(p_buffer[n_x + 1 + n_y * n_width], n_color, n_alpha_0);

					n_lerp_y += n_gradient;
				}
				// draw the line
			} else {
				if(n_end_x0 >= 0 && n_end_x0 < n_width) {
					int n_x = n_end_x0, n_y = int(floor(f_end_y0));
					if(n_y >= 0 && n_y < n_height)
						AlphaBlend(p_buffer[n_x + n_y * n_width], n_color, int((255 - n_alpha_y0) * f_cov_x0));
					if(n_y + 1 >= 0 && n_y + 1 < n_height)
						AlphaBlend(p_buffer[n_x + (n_y + 1) * n_width], n_color, int(n_alpha_y0 * f_cov_x0));
				}
				// handle the first endpoint

				if(n_end_x1 >= 0 && n_end_x1 < n_width) {
					int n_x = n_end_x1, n_y = int(floor(f_end_y1));
					if(n_y >= 0 && n_y < n_height)
						AlphaBlend(p_buffer[n_x + n_y * n_width], n_color, int((255 - n_alpha_y1) * f_cov_x1));
					if(n_y + 1 >= 0 && n_y + 1 < n_height)
						AlphaBlend(p_buffer[n_x + (n_y + 1) * n_width], n_color, int(n_alpha_y1 * f_cov_x1));
				}
				// handle the second endpoint

				for(int n_x = n_end_x0 + 1; n_x < n_end_x1; ++ n_x) {
					int n_y = n_lerp_y >> 8;
					_ASSERTE(n_x >= 0 && n_x < n_width && n_y >= 0 && n_y < n_height);
					int n_alpha_0 = n_lerp_y & 0xff;
					AlphaBlend(p_buffer[n_x + n_y * n_width], n_color, 0xff - n_alpha_0);
					if(n_y + 1 < n_height)
						AlphaBlend(p_buffer[n_x + (n_y + 1) * n_width], n_color, n_alpha_0);

					n_lerp_y += n_gradient;
				}
			}
			// thin lines
		} else {
			float f_dx = fabs(f_x1 - f_x0);
			float f_dy = fabs(f_y1 - f_y0);
			float f_thickness_scale = sqrt(f_dx * f_dx + f_dy * f_dy) / std::max(1.0f, std::max(f_dx, f_dy));
			n_line_width = int(n_line_width * f_thickness_scale + .5f);
			// adjust line thickness based on line angle

			int n_line_extent_top = (n_line_width - 1) / 2;
			int n_line_extent_bottom = (n_line_width - 1) - n_line_extent_top;
			// calculate extent on top and bottom

			float f_end_y0 = f_y0 + f_dxdy * (n_end_x0 - f_x0);
			float f_cov_x0 = ceil(f_x0 + .5f) - (f_x0 + .5f); // how much of line is in the first pixel
			int n_lerp_y = int((f_end_y0 + f_dxdy) * 256);
			float f_end_y1 = f_y1 + f_dxdy * (n_end_x1 - f_x1);
			float f_cov_x1 = 1 - (ceil(f_x1 + .5f) - (f_x1 + .5f)); // how much of line is in the last pixel
			int n_alpha_y0 = 255 - int(255 * (ceil(f_end_y0) - f_end_y0));
			int n_alpha_y1 = int(255 * (f_end_y1 - floor(f_end_y1)));
			// calculate aliassing on the end of the lines
			// note the .5 are important otherwise antialiassing discontinuities occur in the first quadrant

			if(b_steep) {
				if(n_end_x0 >= 0 && n_end_x0 < n_height) {
					int n_y = n_end_x0, n_x = int(floor(f_end_y0));
					if(n_x - n_line_extent_top >= 0 && n_x - n_line_extent_top < n_width)
						AlphaBlend(p_buffer[n_x - n_line_extent_top + n_y * n_width], n_color, int((255 - n_alpha_y0) * f_cov_x0));
					for(int dy = -n_line_extent_top + 1; dy < n_line_extent_bottom; ++ dy) {
						if(n_x + dy >= 0 && n_x + dy < n_width)
							p_buffer[n_x + dy + n_y * n_width] = n_color; //AlphaBlend(p_buffer[n_x + dy + n_y * n_width], n_color, int(255 * f_cov_x0)); // does not give correct results
					}
					if(n_x + n_line_extent_bottom >= 0 && n_x + n_line_extent_bottom < n_width)
						AlphaBlend(p_buffer[n_x + n_line_extent_bottom + n_y * n_width], n_color, int(n_alpha_y0 * f_cov_x0));
				}
				// handle the first endpoint

				if(n_end_x1 >= 0 && n_end_x1 < n_height) {
					int n_y = n_end_x1, n_x = int(floor(f_end_y1));
					if(n_x - n_line_extent_top >= 0 && n_x - n_line_extent_top < n_width)
						AlphaBlend(p_buffer[n_x - n_line_extent_top + n_y * n_width], n_color, int((255 - n_alpha_y1) * f_cov_x1));
					for(int dy = -n_line_extent_top + 1; dy < n_line_extent_bottom; ++ dy) {
						if(n_x + dy >= 0 && n_x + dy < n_width)
							p_buffer[n_x + dy + n_y * n_width] = n_color; //AlphaBlend(p_buffer[n_x + dy + n_y * n_width], n_color, int(255 * f_cov_x1)); // does not give correct results
					}
					if(n_x + n_line_extent_bottom >= 0 && n_x + n_line_extent_bottom < n_width)
						AlphaBlend(p_buffer[n_x + n_line_extent_bottom + n_y * n_width], n_color, int(n_alpha_y1 * f_cov_x1));
				}
				// handle the second endpoint

				for(int n_y = n_end_x0 + 1; n_y < n_end_x1; ++ n_y) {
					int n_x = n_lerp_y >> 8;
					_ASSERTE(n_x >= 0 && n_x < n_width && n_y >= 0 && n_y < n_height);
					int n_alpha_0 = n_lerp_y & 0xff;
					if(n_x >= n_line_extent_top)
						AlphaBlend(p_buffer[n_x - n_line_extent_top + n_y * n_width], n_color, 0xff - n_alpha_0);
					for(int dy = n_line_extent_top - 1; dy > 0; -- dy) {
						if(n_x >= dy)
							p_buffer[n_x - dy + n_y * n_width] = n_color;
					}
					p_buffer[n_x + n_y * n_width] = n_color;
					for(int dy = 1; dy < n_line_extent_bottom; ++ dy) {
						if(n_x + dy < n_width)
							p_buffer[n_x + dy + n_y * n_width] = n_color;
					}
					if(n_x + n_line_extent_bottom < n_width)
						AlphaBlend(p_buffer[n_x + n_line_extent_bottom + n_y * n_width], n_color, n_alpha_0);

					n_lerp_y += n_gradient;
				}
				// draw the line
			} else {
				if(n_end_x0 >= 0 && n_end_x0 < n_width) {
					int n_x = n_end_x0, n_y = int(floor(f_end_y0));
					if(n_y - n_line_extent_top >= 0 && n_y - n_line_extent_top < n_height)
						AlphaBlend(p_buffer[n_x + (n_y - n_line_extent_top) * n_width], n_color, int((255 - n_alpha_y0) * f_cov_x0));
					for(int dy = -n_line_extent_top + 1; dy < n_line_extent_bottom; ++ dy) {
						if(n_y + dy >= 0 && n_y + n_line_extent_bottom < n_height)
							p_buffer[n_x + (n_y + dy) * n_width] = n_color;
					}
					if(n_y + n_line_extent_bottom >= 0 && n_y + n_line_extent_bottom < n_height)
						AlphaBlend(p_buffer[n_x + (n_y + n_line_extent_bottom) * n_width], n_color, int(n_alpha_y0 * f_cov_x0));
				}
				// handle the first endpoint

				if(n_end_x1 >= 0 && n_end_x1 < n_width) {
					int n_x = n_end_x1, n_y = int(floor(f_end_y1));
					if(n_y - n_line_extent_top >= 0 && n_y - n_line_extent_top < n_height)
						AlphaBlend(p_buffer[n_x + (n_y - n_line_extent_top) * n_width], n_color, int((255 - n_alpha_y1) * f_cov_x1));
					for(int dy = -n_line_extent_top + 1; dy < n_line_extent_bottom; ++ dy) {
						if(n_y + dy >= 0 && n_y + n_line_extent_bottom < n_height)
							p_buffer[n_x + (n_y + dy) * n_width] = n_color;
					}
					if(n_y + n_line_extent_bottom >= 0 && n_y + n_line_extent_bottom < n_height)
						AlphaBlend(p_buffer[n_x + (n_y + n_line_extent_bottom) * n_width], n_color, int(n_alpha_y1 * f_cov_x1));
				}
				// handle the second endpoint

				for(int n_x = n_end_x0 + 1; n_x < n_end_x1; ++ n_x) {
					int n_y = n_lerp_y >> 8;
					_ASSERTE(n_x >= 0 && n_x < n_width && n_y >= 0 && n_y < n_height);
					int n_alpha_0 = n_lerp_y & 0xff;
					if(n_y - n_line_extent_top >= 0)
						AlphaBlend(p_buffer[n_x + (n_y - n_line_extent_top) * n_width], n_color, 0xff - n_alpha_0);
					for(int dy = n_line_extent_top - 1; dy > 0; -- dy) {
						if(n_y >= dy)
							p_buffer[n_x + (n_y - dy) * n_width] = n_color;
					}
					p_buffer[n_x + n_y * n_width] = n_color;
					for(int dy = 1; dy < n_line_extent_bottom; ++ dy) {
						if(n_y + dy < n_height)
							p_buffer[n_x + (n_y + dy) * n_width] = n_color;
					}
					if(n_y + n_line_extent_bottom < n_height)
						AlphaBlend(p_buffer[n_x + (n_y + n_line_extent_bottom) * n_width], n_color, n_alpha_0);

					n_lerp_y += n_gradient;
				}
			}
			// thick lines
		}
	}

	/**
	 *	@brief draws a solid-color line
	 *
	 *	@param[in] p_z_buffer is pointer to the buffer containing 1/z values
	 *		(a block of memory with the same dimensions and addressing as this bitmap)
	 *	@param[in] x1 is a coordinate of the first line point
	 *	@param[in] y1 is a coordinate of the first line point
	 *	@param[in] z1 is a coordinate of the first line point
	 *	@param[in] x2 is a coordinate of the second line point
	 *	@param[in] y2 is a coordinate of the second line point
	 *	@param[in] z2 is a coordinate of the second line point
	 *	@param[in] n_color is the line color
	 *	@param[in] n_line_width is the line width, in pixels
	 *
	 *	@note This performs array boundary checking,
	 *		coordinates outside the bitmap are ok.
	 *	@note Depth test behavior is equivalent to that of GL_LESS.
	 *
	 *	@todo Write templated variant of this function to support interpolation of more
	 *		coordinates and custom depth test / shading.
	 *	@todo Implement a better clipping (large offscreen lines
	 *		are rasterized rather slow now).
	 */
	void DrawLine_ZBuffer(float *p_z_buffer, float x1, float y1, float z1,
		float x2, float y2, float z2, uint32_t n_color, int n_line_width = 1)
	{
		float x = floor(x1);
		float y = floor(y1);
		int xlen = int(x1) - int(x2);
		int ylen = int(y1) - int(y2);
		int len = abs(((abs(ylen) > abs(xlen))? ylen : xlen));
		float stepx = float(xlen) / len;
		float stepy = float(ylen) / len;
		float z = 1 / z1;
		float zstep = (1 / z2 - 1 / z1) / len;

		if(n_line_width == 1) {
			for(int i = 0; i <= (int)len; ++ i) {
				if(x >= 0 && x < n_width && y >= 0 && y < n_height && p_z_buffer[int(x) + n_width * int(y)] > z) {
					p_z_buffer[int(x) + n_width * int(y)] = z;
					p_buffer[int(x) + n_width * int(y)] = n_color;
				}
				x -= stepx;
				y -= stepy;
				z += zstep;
			}
		} else {
			int n_start = -n_line_width / 2;
			int n_end = n_start + n_line_width;

			if(abs(xlen) < abs(ylen)) {
				for(int i = 0; i <= (int)len; ++ i) {
					if(y >= 0 && y < n_height) {
						for(int xs = std::max(0, int(x + n_start)), xe = std::min(n_width, int(x + n_end)); xs < xe; ++ xs) {
							if(p_z_buffer[xs + n_width * int(y)] > z) {
								p_z_buffer[xs + n_width * int(y)] = z;
								p_buffer[xs + n_width * int(y)] = n_color;
							}
						}
					}
					x -= stepx;
					y -= stepy;
					z += zstep;
				}
			} else {
				for(int i = 0; i <= (int)len; ++ i) {
					if(x >= 0 && x < n_width) {
						for(int ys = std::max(0, int(y + n_start)), ye = std::min(n_height, int(y + n_end)); ys < ye; ++ ys) {
							if(p_z_buffer[int(x) + n_width * ys] > z) {
								p_z_buffer[int(x) + n_width * ys] = z;
								p_buffer[int(x) + n_width * ys] = n_color;
							}
						}
					}
					x -= stepx;
					y -= stepy;
					z += zstep;
				}
			}
		}
		// quick and dirty
	}

	/**
	 *	@brief draws a solid-color triangle
	 *
	 *	@param[in] p_z_buffer is pointer to the buffer containing 1/z values
	 *		(a block of memory with the same dimensions and addressing as this bitmap)
	 *	@param[in] p_vertex is a list of the triangle vertices
	 *	@param[in] n_color is the fill color
	 *
	 *	@note This performs clipping, coordinates outside the bitmap are ok.
	 *	@note Depth test behavior is equivalent to that of GL_LESS.
	 *
	 *	@todo Write templated variant of this function to support interpolation of more
	 *		coordinates and custom depth test / shading.
	 */
	void DrawTriangle_ZBuffer(float *p_z_buffer, const float p_vertex[3][3], uint32_t n_color)
	{
		int n_min_y = int(floor(std::min(p_vertex[0][1], std::min(p_vertex[1][1], p_vertex[2][1]))));
		int n_max_y = int(floor(std::max(p_vertex[0][1], std::max(p_vertex[1][1], p_vertex[2][1]))));
		// find the top / bottom y

		if(n_min_y < 0)
			n_min_y = 0;
		if(n_max_y >= n_height)
			n_max_y = n_height - 1;
		// make sure we don't compute pixels that are not displayed in the end

		p_z_buffer += n_width * n_min_y;
		uint32_t *p_color_buffer = p_buffer + n_width * n_min_y;
		for(int y = n_min_y; y <= n_max_y; ++ y, p_z_buffer += n_width, p_color_buffer += n_width) {
			int n_point_num = 0;
			float p_segment[2][3];
			for(int i = 0, j = 2; i < 3; j = i ++) {
				const float *p_a = p_vertex[i];
				const float *p_b = p_vertex[j];

				if((y >= p_a[1] && y <= p_b[1]) ||
				   (y >= p_b[1] && y <= p_a[1])) {
					float t = (y - p_a[1]) / (p_b[1] - p_a[1]);
					// find edge-scanline intersection

					if(n_point_num < 2) {
						for(int n = 0; n < 3; ++ n)
							p_segment[n_point_num][n] = p_a[n] + t * (p_b[n] - p_a[n]);
						++ n_point_num;
					} else {
						if(p_segment[0][0] > p_segment[1][0]) {
							for(int n = 0; n < 3; ++ n)
								std::swap(p_segment[0][n], p_segment[1][n]);
						}
						// make sure the first is left and the second is right

						float p_isect[3];
						for(int n = 0; n < 3; ++ n)
							p_isect[n] = p_a[n] + t * (p_b[n] - p_a[n]);
						// calculate the new intersection

						if(p_isect[0] < p_segment[0][0])
							memcpy(p_segment[0], p_isect, 3 * sizeof(float)); // new is left to 0
						else if(p_isect[0] > p_segment[1][0])
							memcpy(p_segment[1], p_isect, 3 * sizeof(float)); // new is right to 1
						// in case the third intersection was found on the same scanline,
						// replace the one lying in the segment of the other two
					}
					// calculate intersection position, add it to the list
				}
			}
			// find intersection of the triangle and the scanline

			if(n_point_num != 2)
				continue;
			// bad intersections

			if(p_segment[0][0] > p_segment[1][0]) {
				for(int n = 0; n < 3; ++ n)
					std::swap(p_segment[0][n], p_segment[1][n]);
			}
			// make sure the first is left and the second is right

			if(int(p_segment[1][0]) < 0 || p_segment[0][0] >= n_width)
				continue;
			// it's too left, or too right

			p_segment[0][2] = 1 / p_segment[0][2];
			p_segment[1][2] = 1 / p_segment[1][2];
			// convert z to 1/z to make it's interpolation linear

			float p_delta[3];
			for(int n = 0; n < 3; ++ n)
				p_delta[n] = p_segment[1][n] - p_segment[0][n];
			float f_len = p_delta[0];
			for(int m = 0; m < 3; ++ m)
				p_delta[m] /= f_len;
			// calculate delta coordinates per x-step

			int l = std::max(0, int(floor(p_segment[0][0])));
			if(p_segment[0][0] != l) {
				float f_offset = l - p_segment[0][0];
				for(int n = 0; n < 3; ++ n)
					p_segment[0][n] += p_delta[n] * f_offset;
			}
			// fixup left point if offscreen

			int r = std::min(n_width - 1, int(floor(p_segment[1][0])));
			for(; l <= r; ++ l) {
				if(p_segment[0][2] < p_z_buffer[l]) {
					p_z_buffer[l] = p_segment[0][2];
					p_color_buffer[l] = n_color;
				}

				for(int n = 0; n < 3; ++ n)
					p_segment[0][n] += p_delta[n];
			}
			// rasterize the segment
		}
		// rasterize the triangle
	}
};

#endif // __BITMAP_STRUCTURE_INCLUDED
