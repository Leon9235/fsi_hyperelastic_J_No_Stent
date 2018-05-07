/**
 * @file BioChem.h
 * @brief 
 * @author Shihua Gong
 * @version 
 * @date 30-04-2016
 */

#ifndef __Common_H__
#define __Common_H__ 

enum FacetType : std::size_t { GAMMA_FI = 202, GAMMA_FO = 203,
							   GAMMA_F2I = 211, GAMMA_F2O = 212,
							   GAMMA_VI = 215, GAMMA_VO = 216, GAMMA_VOUT = 217,
							   GAMMA_SI = 205, GAMMA_SO = 206,
							   GAMMA_IFS = 201, GAMMA_IFV = 219, GAMMA_IF2V = 213};

//LZ GAMMA_FI inlet flow of BLOOD
//LZ GAMMA_FO outlet flow  of BLOOD
//LZ GAMMA_F2I inlet flow of BLOOD2
//LZ GAMMA_F2O outlet flow  of BLOOD2
//LZ GAMMA_VI inlet of VESSEL
//LZ GAMMA_VO outlet of VESSEL
//LZ GAMMA_VOUT outer wall of VESSEL
//LZ GAMMA_SI inlet of STENT
//LZ GAMMA_SO outlet of STENT
//LZ GAMMA_IFS interface of BLOOD and STENT
//LZ GAMMA_IFV interface of BLOOD and VESSEL
//LZ GAMMA_IF2V interface of BLOOD2 and VESSEL


enum DomainType : std::size_t {BLOOD = 204, BLOOD2 = 214, VESSEL = 218, STENT = 210}; 

//LZ BLOOD blood inside the stent + the other parts 
//LZ BLOOD2 blood surrounding the stent 
//LZ VESSEL 
//LZ STENT

#endif // __Common_H__

/** 
 * end of file 
 *
 */
