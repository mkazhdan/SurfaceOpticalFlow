/*
Copyright (c) 2010, Michael Kazhdan
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this list of
conditions and the following disclaimer. Redistributions in binary form must reproduce
the above copyright notice, this list of conditions and the following disclaimer
in the documentation and/or other materials provided with the distribution. 

Neither the name of the Johns Hopkins University nor the names of its contributors
may be used to endorse or promote products derived from this software without specific
prior written permission. 

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO THE IMPLIED WARRANTIES 
OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
TO, PROCUREMENT OF SUBSTITUTE  GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.
*/

#include <stdio.h>
#include <stdlib.h>
#include <setjmp.h>

#ifdef _WIN32
#include <windows.h>
#include "JPEG/jpeglib.h"
#include "JPEG/jerror.h"
#include "JPEG/jmorecfg.h"
#else // !_WIN32
#include <jpeglib.h>
#include <jerror.h>
#include <jmorecfg.h>
#endif // _WIN32

struct my_error_mgr
{
	struct jpeg_error_mgr pub;    // "public" fields
	jmp_buf setjmp_buffer;        // for return to caller
};
typedef struct my_error_mgr * my_error_ptr;


METHODDEF( void )
my_error_exit (j_common_ptr cinfo)
{
	// cinfo->err really points to a my_error_mgr struct, so coerce pointer
	my_error_ptr myerr = (my_error_ptr) cinfo->err;

	// Always display the message.
	// We could postpone this until after returning, if we chose.
	(*cinfo->err->output_message) (cinfo);

	// Return control to the setjmp point
	longjmp(myerr->setjmp_buffer, 1);
}

void JPEGWriteColor( const char* fileName , unsigned char* pixels , int width , int height , int quality )
{
	FILE* fp;
	struct jpeg_compress_struct cInfo;
	struct my_error_mgr jErr;

	fp = fopen( fileName , "wb" );
	if( !fp ) fprintf( stderr , "[ERROR] Failed to open: %s\n" , fileName ) , exit(0);

	cInfo.err = jpeg_std_error( &jErr.pub );
	jpeg_create_compress( &cInfo );

	jpeg_stdio_dest( &cInfo , fp );

	cInfo.image_width = width;			/* image width and height, in pixels */
	cInfo.image_height = height;
	cInfo.input_components = 3;			/* # of color components per pixel */
	cInfo.in_color_space = JCS_RGB;		/* colorspace of input image */

	jpeg_set_defaults( &cInfo );
	jpeg_set_quality( &cInfo , quality , TRUE );

	jpeg_start_compress( &cInfo , TRUE );

	for( int j=0 ; j<height ; j++ )
	{
		JSAMPROW row_pointer[1];
		row_pointer[0] = pixels + j * width * 3;
		(void) jpeg_write_scanlines( &cInfo , row_pointer , 1 );
	}
	jpeg_finish_compress( &cInfo );
	jpeg_destroy_compress( &cInfo );
	fclose( fp );
}

unsigned char* JPEGReadColor( const char* fileName , int& width , int& height )
{
	FILE* fp;
	struct jpeg_decompress_struct cInfo;
	struct my_error_mgr jErr;

	fp = fopen( fileName , "rb" );
	if( !fp ) fprintf( stderr , "[ERROR] Failed to open: %s\n" , fileName ) , exit(0);

	cInfo.err = jpeg_std_error( &jErr.pub );
	jErr.pub.error_exit = my_error_exit;
	if( setjmp( jErr.setjmp_buffer ) )
	{
		jpeg_destroy_decompress( &cInfo );
		fprintf( stderr , "[ERROR] JPEG error occured\n" );
		exit( 0 );
	}

	jpeg_create_decompress( &cInfo );
	jpeg_stdio_src( &cInfo , fp );

	(void) jpeg_read_header( &cInfo , TRUE );
	(void) jpeg_start_decompress( &cInfo );

	if( cInfo.output_components!=3 )
	{
		fprintf( stderr , "[ERROR] Only 3 components per pixel supported: %d != 3\n" , cInfo.output_components );
		exit( 0 );
	}
	width  = cInfo.output_width;
	height = cInfo.output_height;

	unsigned char* pixels = (unsigned char*) malloc( sizeof( unsigned char ) * width * height * 3 );
	if( !pixels )
	{
		fprintf( stderr , "[ERROR] Failed to allocate pixels: %d x %d\n" , width , height );
		exit( 0 );
	}

	for( int j=0 ; j<height ; j++ )
	{
		JSAMPROW row_pointers[1];
		row_pointers[0] = pixels + j * width * 3;
		jpeg_read_scanlines( &cInfo , row_pointers, 1 );
	}
	(void) jpeg_finish_decompress( &cInfo );
	jpeg_destroy_decompress( &cInfo );
	fclose( fp );
	return pixels;
}