
typedef   signed char      ssp_s8;
typedef unsigned char      ssp_u8;
 
typedef   signed short     ssp_s16;
typedef unsigned short     ssp_u16;
 
typedef   signed int       ssp_s32;
typedef unsigned int       ssp_u32;
 
typedef float              ssp_f32;
typedef double             ssp_f64;
 
typedef   signed long long ssp_s64;
typedef unsigned long long ssp_u64;
 
typedef union 
{
  __m128  f;
  __m128d d;
  __m128i i;
  __m64       m64[ 2];
  ssp_u64 u64[ 2];  
  ssp_s64 s64[ 2];
  ssp_f64 f64[ 2]; 
  ssp_u32 u32[ 4];
  ssp_s32 s32[ 4];    
  ssp_f32 f32[ 4]; 
  ssp_u16 u16[ 8];  
  ssp_s16 s16[ 8];    
  ssp_u8  u8 [16];
  ssp_s8  s8 [16];    
} ssp_m128;

inline __m128i _mm_abs_epi32 	( __m128i  a  )   	
{
  __m128i mask = _mm_cmplt_epi32( a, _mm_setzero_si128() ); 
  a    = _mm_xor_si128 ( a, mask );                        
  mask = _mm_srli_epi32( mask, 31 );                       
  a = _mm_add_epi32( a, mask );                           
  return a;
}


inline void ssp_convert_odd_even_ps_SSE2  	(  __m128 *   	 a,
						   __m128 *  	b	 
						   ) 			
{
  __m128 c, d;  
  c = _mm_shuffle_ps( *a, *b, _MM_SHUFFLE(3,1,3,1) );
  d = _mm_shuffle_ps( *a, *b, _MM_SHUFFLE(2,0,2,0) );
  *a = c;
  *b = d;     
}

inline  void ssp_convert_odd_even_epi32_SSE2  	(  	__m128i *   	 a,
							__m128i *  	b	 
							) 			
{
  ssp_m128 A,B;
  A.i = *a;
  B.i = *b;  
   
  ssp_convert_odd_even_ps_SSE2( &A.f, &B.f );
   
  *a = A.i;
  *b = B.i;       
}

inline __m128i _mm_hadd_epi32  	(__m128i   	 a,
				 __m128i  	b	 
				 ) 			
{
  ssp_convert_odd_even_epi32_SSE2( &a, &b );
  a = _mm_add_epi32( a, b );
  return a; 
}



