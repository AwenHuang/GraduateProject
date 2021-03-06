﻿; Listing generated by Microsoft (R) Optimizing Compiler Version 16.00.30319.01 

	TITLE	D:\專題\自建專案\自建專案\src\BulletSoftBody\btSoftSoftCollisionAlgorithm.cpp
	.686P
	.XMM
	include listing.inc
	.model	flat

INCLUDELIB LIBCMTD
INCLUDELIB OLDNAMES

_DATA	SEGMENT
_btNanMask DD	07f800001H
_btInfinityMask DD 07f800000H
_DATA	ENDS
PUBLIC	??_H@YGXPAXIHP6EPAX0@Z@Z			; `vector constructor iterator'
EXTRN	__RTC_CheckEsp:PROC
EXTRN	__RTC_Shutdown:PROC
EXTRN	__RTC_InitBase:PROC
;	COMDAT rtc$TMZ
rtc$TMZ	SEGMENT
__RTC_Shutdown.rtc$TMZ DD FLAT:__RTC_Shutdown
rtc$TMZ	ENDS
;	COMDAT rtc$IMZ
rtc$IMZ	SEGMENT
__RTC_InitBase.rtc$IMZ DD FLAT:__RTC_InitBase
; Function compile flags: /Odtp /RTCsu
rtc$IMZ	ENDS
;	COMDAT ??_H@YGXPAXIHP6EPAX0@Z@Z
_TEXT	SEGMENT
___t$ = 8						; size = 4
___s$ = 12						; size = 4
___n$ = 16						; size = 4
___f$ = 20						; size = 4
??_H@YGXPAXIHP6EPAX0@Z@Z PROC				; `vector constructor iterator', COMDAT
	push	ebp
	mov	ebp, esp
	push	esi
$LN2@vector:
	mov	eax, DWORD PTR ___n$[ebp]
	sub	eax, 1
	mov	DWORD PTR ___n$[ebp], eax
	js	SHORT $LN3@vector
	mov	esi, esp
	mov	ecx, DWORD PTR ___t$[ebp]
	call	DWORD PTR ___f$[ebp]
	cmp	esi, esp
	call	__RTC_CheckEsp
	mov	ecx, DWORD PTR ___t$[ebp]
	add	ecx, DWORD PTR ___s$[ebp]
	mov	DWORD PTR ___t$[ebp], ecx
	jmp	SHORT $LN2@vector
$LN3@vector:
	pop	esi
	cmp	ebp, esp
	call	__RTC_CheckEsp
	pop	ebp
	ret	16					; 00000010H
??_H@YGXPAXIHP6EPAX0@Z@Z ENDP				; `vector constructor iterator'
_TEXT	ENDS
PUBLIC	??_7btSoftSoftCollisionAlgorithm@@6B@		; btSoftSoftCollisionAlgorithm::`vftable'
PUBLIC	??0btSoftSoftCollisionAlgorithm@@QAE@PAVbtPersistentManifold@@ABUbtCollisionAlgorithmConstructionInfo@@PBUbtCollisionObjectWrapper@@2@Z ; btSoftSoftCollisionAlgorithm::btSoftSoftCollisionAlgorithm
PUBLIC	?processCollision@btSoftSoftCollisionAlgorithm@@UAEXPBUbtCollisionObjectWrapper@@0ABUbtDispatcherInfo@@PAVbtManifoldResult@@@Z ; btSoftSoftCollisionAlgorithm::processCollision
PUBLIC	?calculateTimeOfImpact@btSoftSoftCollisionAlgorithm@@UAEMPAVbtCollisionObject@@0ABUbtDispatcherInfo@@PAVbtManifoldResult@@@Z ; btSoftSoftCollisionAlgorithm::calculateTimeOfImpact
PUBLIC	?getAllContactManifolds@btSoftSoftCollisionAlgorithm@@UAEXAAV?$btAlignedObjectArray@PAVbtPersistentManifold@@@@@Z ; btSoftSoftCollisionAlgorithm::getAllContactManifolds
EXTRN	??0btCollisionAlgorithm@@QAE@ABUbtCollisionAlgorithmConstructionInfo@@@Z:PROC ; btCollisionAlgorithm::btCollisionAlgorithm
EXTRN	??_EbtSoftSoftCollisionAlgorithm@@UAEPAXI@Z:PROC ; btSoftSoftCollisionAlgorithm::`vector deleting destructor'
;	COMDAT ??_7btSoftSoftCollisionAlgorithm@@6B@
; File d:\專題\自建專案\自建專案\src\bulletsoftbody\btsoftsoftcollisionalgorithm.cpp
CONST	SEGMENT
??_7btSoftSoftCollisionAlgorithm@@6B@ DD FLAT:??_EbtSoftSoftCollisionAlgorithm@@UAEPAXI@Z ; btSoftSoftCollisionAlgorithm::`vftable'
	DD	FLAT:?processCollision@btSoftSoftCollisionAlgorithm@@UAEXPBUbtCollisionObjectWrapper@@0ABUbtDispatcherInfo@@PAVbtManifoldResult@@@Z
	DD	FLAT:?calculateTimeOfImpact@btSoftSoftCollisionAlgorithm@@UAEMPAVbtCollisionObject@@0ABUbtDispatcherInfo@@PAVbtManifoldResult@@@Z
	DD	FLAT:?getAllContactManifolds@btSoftSoftCollisionAlgorithm@@UAEXAAV?$btAlignedObjectArray@PAVbtPersistentManifold@@@@@Z
; Function compile flags: /Odtp /RTCsu
CONST	ENDS
;	COMDAT ??0btSoftSoftCollisionAlgorithm@@QAE@PAVbtPersistentManifold@@ABUbtCollisionAlgorithmConstructionInfo@@PBUbtCollisionObjectWrapper@@2@Z
_TEXT	SEGMENT
_this$ = -4						; size = 4
___formal$ = 8						; size = 4
_ci$ = 12						; size = 4
___formal$ = 16						; size = 4
___formal$ = 20						; size = 4
??0btSoftSoftCollisionAlgorithm@@QAE@PAVbtPersistentManifold@@ABUbtCollisionAlgorithmConstructionInfo@@PBUbtCollisionObjectWrapper@@2@Z PROC ; btSoftSoftCollisionAlgorithm::btSoftSoftCollisionAlgorithm, COMDAT
; _this$ = ecx
; Line 30
	push	ebp
	mov	ebp, esp
	push	ecx
	mov	DWORD PTR [ebp-4], -858993460		; ccccccccH
	mov	DWORD PTR _this$[ebp], ecx
	mov	eax, DWORD PTR _ci$[ebp]
	push	eax
	mov	ecx, DWORD PTR _this$[ebp]
	call	??0btCollisionAlgorithm@@QAE@ABUbtCollisionAlgorithmConstructionInfo@@@Z ; btCollisionAlgorithm::btCollisionAlgorithm
	mov	ecx, DWORD PTR _this$[ebp]
	mov	DWORD PTR [ecx], OFFSET ??_7btSoftSoftCollisionAlgorithm@@6B@
; Line 31
	mov	eax, DWORD PTR _this$[ebp]
	add	esp, 4
	cmp	ebp, esp
	call	__RTC_CheckEsp
	mov	esp, ebp
	pop	ebp
	ret	16					; 00000010H
??0btSoftSoftCollisionAlgorithm@@QAE@PAVbtPersistentManifold@@ABUbtCollisionAlgorithmConstructionInfo@@PBUbtCollisionObjectWrapper@@2@Z ENDP ; btSoftSoftCollisionAlgorithm::btSoftSoftCollisionAlgorithm
_TEXT	ENDS
PUBLIC	?push_back@?$btAlignedObjectArray@PAVbtPersistentManifold@@@@QAEXABQAVbtPersistentManifold@@@Z ; btAlignedObjectArray<btPersistentManifold *>::push_back
; Function compile flags: /Odtp /RTCsu
;	COMDAT ?getAllContactManifolds@btSoftSoftCollisionAlgorithm@@UAEXAAV?$btAlignedObjectArray@PAVbtPersistentManifold@@@@@Z
_TEXT	SEGMENT
_this$ = -4						; size = 4
_manifoldArray$ = 8					; size = 4
?getAllContactManifolds@btSoftSoftCollisionAlgorithm@@UAEXAAV?$btAlignedObjectArray@PAVbtPersistentManifold@@@@@Z PROC ; btSoftSoftCollisionAlgorithm::getAllContactManifolds, COMDAT
; _this$ = ecx
; File d:\專題\自建專案\自建專案\src\bulletsoftbody\btsoftsoftcollisionalgorithm.h
; Line 46
	push	ebp
	mov	ebp, esp
	push	ecx
	mov	DWORD PTR [ebp-4], -858993460		; ccccccccH
	mov	DWORD PTR _this$[ebp], ecx
; Line 47
	mov	eax, DWORD PTR _this$[ebp]
	cmp	DWORD PTR [eax+12], 0
	je	SHORT $LN2@getAllCont
	mov	ecx, DWORD PTR _this$[ebp]
	movzx	edx, BYTE PTR [ecx+8]
	test	edx, edx
	je	SHORT $LN2@getAllCont
; Line 48
	mov	eax, DWORD PTR _this$[ebp]
	add	eax, 12					; 0000000cH
	push	eax
	mov	ecx, DWORD PTR _manifoldArray$[ebp]
	call	?push_back@?$btAlignedObjectArray@PAVbtPersistentManifold@@@@QAEXABQAVbtPersistentManifold@@@Z ; btAlignedObjectArray<btPersistentManifold *>::push_back
$LN2@getAllCont:
; Line 49
	add	esp, 4
	cmp	ebp, esp
	call	__RTC_CheckEsp
	mov	esp, ebp
	pop	ebp
	ret	4
?getAllContactManifolds@btSoftSoftCollisionAlgorithm@@UAEXAAV?$btAlignedObjectArray@PAVbtPersistentManifold@@@@@Z ENDP ; btSoftSoftCollisionAlgorithm::getAllContactManifolds
_TEXT	ENDS
PUBLIC	??1btSoftSoftCollisionAlgorithm@@UAE@XZ		; btSoftSoftCollisionAlgorithm::~btSoftSoftCollisionAlgorithm
EXTRN	??3@YAXPAX@Z:PROC				; operator delete
; Function compile flags: /Odtp /RTCsu
;	COMDAT ??_GbtSoftSoftCollisionAlgorithm@@UAEPAXI@Z
_TEXT	SEGMENT
_this$ = -4						; size = 4
___flags$ = 8						; size = 4
??_GbtSoftSoftCollisionAlgorithm@@UAEPAXI@Z PROC	; btSoftSoftCollisionAlgorithm::`scalar deleting destructor', COMDAT
; _this$ = ecx
	push	ebp
	mov	ebp, esp
	push	ecx
	mov	DWORD PTR [ebp-4], -858993460		; ccccccccH
	mov	DWORD PTR _this$[ebp], ecx
	mov	ecx, DWORD PTR _this$[ebp]
	call	??1btSoftSoftCollisionAlgorithm@@UAE@XZ	; btSoftSoftCollisionAlgorithm::~btSoftSoftCollisionAlgorithm
	mov	eax, DWORD PTR ___flags$[ebp]
	and	eax, 1
	je	SHORT $LN1@scalar
	mov	ecx, DWORD PTR _this$[ebp]
	push	ecx
	call	??3@YAXPAX@Z				; operator delete
	add	esp, 4
$LN1@scalar:
	mov	eax, DWORD PTR _this$[ebp]
	add	esp, 4
	cmp	ebp, esp
	call	__RTC_CheckEsp
	mov	esp, ebp
	pop	ebp
	ret	4
??_GbtSoftSoftCollisionAlgorithm@@UAEPAXI@Z ENDP	; btSoftSoftCollisionAlgorithm::`scalar deleting destructor'
_TEXT	ENDS
PUBLIC	??1btCollisionAlgorithm@@UAE@XZ			; btCollisionAlgorithm::~btCollisionAlgorithm
; Function compile flags: /Odtp /RTCsu
;	COMDAT ??1btSoftSoftCollisionAlgorithm@@UAE@XZ
_TEXT	SEGMENT
_this$ = -4						; size = 4
??1btSoftSoftCollisionAlgorithm@@UAE@XZ PROC		; btSoftSoftCollisionAlgorithm::~btSoftSoftCollisionAlgorithm, COMDAT
; _this$ = ecx
; File d:\專題\自建專案\自建專案\src\bulletsoftbody\btsoftsoftcollisionalgorithm.cpp
; Line 34
	push	ebp
	mov	ebp, esp
	push	ecx
	mov	DWORD PTR [ebp-4], -858993460		; ccccccccH
	mov	DWORD PTR _this$[ebp], ecx
	mov	eax, DWORD PTR _this$[ebp]
	mov	DWORD PTR [eax], OFFSET ??_7btSoftSoftCollisionAlgorithm@@6B@
; Line 35
	mov	ecx, DWORD PTR _this$[ebp]
	call	??1btCollisionAlgorithm@@UAE@XZ		; btCollisionAlgorithm::~btCollisionAlgorithm
	add	esp, 4
	cmp	ebp, esp
	call	__RTC_CheckEsp
	mov	esp, ebp
	pop	ebp
	ret	0
??1btSoftSoftCollisionAlgorithm@@UAE@XZ ENDP		; btSoftSoftCollisionAlgorithm::~btSoftSoftCollisionAlgorithm
_TEXT	ENDS
PUBLIC	??_7btCollisionAlgorithm@@6B@			; btCollisionAlgorithm::`vftable'
EXTRN	??_EbtCollisionAlgorithm@@UAEPAXI@Z:PROC	; btCollisionAlgorithm::`vector deleting destructor'
EXTRN	__purecall:PROC
;	COMDAT ??_7btCollisionAlgorithm@@6B@
; File d:\專題\自建專案\自建專案\src\bulletcollision\broadphasecollision\btcollisionalgorithm.h
CONST	SEGMENT
??_7btCollisionAlgorithm@@6B@ DD FLAT:??_EbtCollisionAlgorithm@@UAEPAXI@Z ; btCollisionAlgorithm::`vftable'
	DD	FLAT:__purecall
	DD	FLAT:__purecall
	DD	FLAT:__purecall
; Function compile flags: /Odtp /RTCsu
CONST	ENDS
;	COMDAT ??1btCollisionAlgorithm@@UAE@XZ
_TEXT	SEGMENT
_this$ = -4						; size = 4
??1btCollisionAlgorithm@@UAE@XZ PROC			; btCollisionAlgorithm::~btCollisionAlgorithm, COMDAT
; _this$ = ecx
; Line 71
	push	ebp
	mov	ebp, esp
	push	ecx
	mov	DWORD PTR [ebp-4], -858993460		; ccccccccH
	mov	DWORD PTR _this$[ebp], ecx
	mov	eax, DWORD PTR _this$[ebp]
	mov	DWORD PTR [eax], OFFSET ??_7btCollisionAlgorithm@@6B@
	mov	esp, ebp
	pop	ebp
	ret	0
??1btCollisionAlgorithm@@UAE@XZ ENDP			; btCollisionAlgorithm::~btCollisionAlgorithm
; Function compile flags: /Odtp /RTCsu
_TEXT	ENDS
;	COMDAT ??_GbtCollisionAlgorithm@@UAEPAXI@Z
_TEXT	SEGMENT
_this$ = -4						; size = 4
___flags$ = 8						; size = 4
??_GbtCollisionAlgorithm@@UAEPAXI@Z PROC		; btCollisionAlgorithm::`scalar deleting destructor', COMDAT
; _this$ = ecx
	push	ebp
	mov	ebp, esp
	push	ecx
	mov	DWORD PTR [ebp-4], -858993460		; ccccccccH
	mov	DWORD PTR _this$[ebp], ecx
	mov	ecx, DWORD PTR _this$[ebp]
	call	??1btCollisionAlgorithm@@UAE@XZ		; btCollisionAlgorithm::~btCollisionAlgorithm
	mov	eax, DWORD PTR ___flags$[ebp]
	and	eax, 1
	je	SHORT $LN1@scalar@2
	mov	ecx, DWORD PTR _this$[ebp]
	push	ecx
	call	??3@YAXPAX@Z				; operator delete
	add	esp, 4
$LN1@scalar@2:
	mov	eax, DWORD PTR _this$[ebp]
	add	esp, 4
	cmp	ebp, esp
	call	__RTC_CheckEsp
	mov	esp, ebp
	pop	ebp
	ret	4
??_GbtCollisionAlgorithm@@UAEPAXI@Z ENDP		; btCollisionAlgorithm::`scalar deleting destructor'
_TEXT	ENDS
PUBLIC	?getSoftBodySolver@btSoftBody@@QAEPAVbtSoftBodySolver@@XZ ; btSoftBody::getSoftBodySolver
PUBLIC	?getCollisionObject@btCollisionObjectWrapper@@QBEPBVbtCollisionObject@@XZ ; btCollisionObjectWrapper::getCollisionObject
; Function compile flags: /Odtp /RTCsu
;	COMDAT ?processCollision@btSoftSoftCollisionAlgorithm@@UAEXPBUbtCollisionObjectWrapper@@0ABUbtDispatcherInfo@@PAVbtManifoldResult@@@Z
_TEXT	SEGMENT
tv79 = -16						; size = 4
_soft1$ = -12						; size = 4
_soft0$ = -8						; size = 4
_this$ = -4						; size = 4
_body0Wrap$ = 8						; size = 4
_body1Wrap$ = 12					; size = 4
___formal$ = 16						; size = 4
___formal$ = 20						; size = 4
?processCollision@btSoftSoftCollisionAlgorithm@@UAEXPBUbtCollisionObjectWrapper@@0ABUbtDispatcherInfo@@PAVbtManifoldResult@@@Z PROC ; btSoftSoftCollisionAlgorithm::processCollision, COMDAT
; _this$ = ecx
; File d:\專題\自建專案\自建專案\src\bulletsoftbody\btsoftsoftcollisionalgorithm.cpp
; Line 38
	push	ebp
	mov	ebp, esp
	sub	esp, 16					; 00000010H
	push	esi
	mov	eax, -858993460				; ccccccccH
	mov	DWORD PTR [ebp-16], eax
	mov	DWORD PTR [ebp-12], eax
	mov	DWORD PTR [ebp-8], eax
	mov	DWORD PTR [ebp-4], eax
	mov	DWORD PTR _this$[ebp], ecx
; Line 39
	mov	ecx, DWORD PTR _body0Wrap$[ebp]
	call	?getCollisionObject@btCollisionObjectWrapper@@QBEPBVbtCollisionObject@@XZ ; btCollisionObjectWrapper::getCollisionObject
	mov	DWORD PTR _soft0$[ebp], eax
; Line 40
	mov	ecx, DWORD PTR _body1Wrap$[ebp]
	call	?getCollisionObject@btCollisionObjectWrapper@@QBEPBVbtCollisionObject@@XZ ; btCollisionObjectWrapper::getCollisionObject
	mov	DWORD PTR _soft1$[ebp], eax
; Line 41
	mov	ecx, DWORD PTR _soft0$[ebp]
	call	?getSoftBodySolver@btSoftBody@@QAEPAVbtSoftBodySolver@@XZ ; btSoftBody::getSoftBodySolver
	mov	DWORD PTR tv79[ebp], eax
	mov	esi, esp
	mov	eax, DWORD PTR _soft1$[ebp]
	push	eax
	mov	ecx, DWORD PTR _soft0$[ebp]
	push	ecx
	mov	edx, DWORD PTR tv79[ebp]
	mov	eax, DWORD PTR [edx]
	mov	ecx, DWORD PTR tv79[ebp]
	mov	edx, DWORD PTR [eax+32]
	call	edx
	cmp	esi, esp
	call	__RTC_CheckEsp
; Line 42
	pop	esi
	add	esp, 16					; 00000010H
	cmp	ebp, esp
	call	__RTC_CheckEsp
	mov	esp, ebp
	pop	ebp
	ret	16					; 00000010H
?processCollision@btSoftSoftCollisionAlgorithm@@UAEXPBUbtCollisionObjectWrapper@@0ABUbtDispatcherInfo@@PAVbtManifoldResult@@@Z ENDP ; btSoftSoftCollisionAlgorithm::processCollision
; Function compile flags: /Odtp /RTCsu
_TEXT	ENDS
;	COMDAT ?getCollisionObject@btCollisionObjectWrapper@@QBEPBVbtCollisionObject@@XZ
_TEXT	SEGMENT
_this$ = -4						; size = 4
?getCollisionObject@btCollisionObjectWrapper@@QBEPBVbtCollisionObject@@XZ PROC ; btCollisionObjectWrapper::getCollisionObject, COMDAT
; _this$ = ecx
; File d:\專題\自建專案\自建專案\src\bulletcollision\collisiondispatch\btcollisionobjectwrapper.h
; Line 39
	push	ebp
	mov	ebp, esp
	push	ecx
	mov	DWORD PTR [ebp-4], -858993460		; ccccccccH
	mov	DWORD PTR _this$[ebp], ecx
	mov	eax, DWORD PTR _this$[ebp]
	mov	eax, DWORD PTR [eax+8]
	mov	esp, ebp
	pop	ebp
	ret	0
?getCollisionObject@btCollisionObjectWrapper@@QBEPBVbtCollisionObject@@XZ ENDP ; btCollisionObjectWrapper::getCollisionObject
; Function compile flags: /Odtp /RTCsu
_TEXT	ENDS
;	COMDAT ?getSoftBodySolver@btSoftBody@@QAEPAVbtSoftBodySolver@@XZ
_TEXT	SEGMENT
_this$ = -4						; size = 4
?getSoftBodySolver@btSoftBody@@QAEPAVbtSoftBodySolver@@XZ PROC ; btSoftBody::getSoftBodySolver, COMDAT
; _this$ = ecx
; File d:\專題\自建專案\自建專案\src\bulletsoftbody\btsoftbody.h
; Line 920
	push	ebp
	mov	ebp, esp
	push	ecx
	mov	DWORD PTR [ebp-4], -858993460		; ccccccccH
	mov	DWORD PTR _this$[ebp], ecx
; Line 921
	mov	eax, DWORD PTR _this$[ebp]
	mov	eax, DWORD PTR [eax+308]
; Line 922
	mov	esp, ebp
	pop	ebp
	ret	0
?getSoftBodySolver@btSoftBody@@QAEPAVbtSoftBodySolver@@XZ ENDP ; btSoftBody::getSoftBodySolver
_TEXT	ENDS
PUBLIC	__real@3f800000
EXTRN	__fltused:DWORD
;	COMDAT __real@3f800000
; File d:\專題\自建專案\自建專案\src\bulletsoftbody\btsoftsoftcollisionalgorithm.cpp
CONST	SEGMENT
__real@3f800000 DD 03f800000r			; 1
; Function compile flags: /Odtp /RTCsu
CONST	ENDS
;	COMDAT ?calculateTimeOfImpact@btSoftSoftCollisionAlgorithm@@UAEMPAVbtCollisionObject@@0ABUbtDispatcherInfo@@PAVbtManifoldResult@@@Z
_TEXT	SEGMENT
_this$ = -4						; size = 4
___formal$ = 8						; size = 4
___formal$ = 12						; size = 4
___formal$ = 16						; size = 4
___formal$ = 20						; size = 4
?calculateTimeOfImpact@btSoftSoftCollisionAlgorithm@@UAEMPAVbtCollisionObject@@0ABUbtDispatcherInfo@@PAVbtManifoldResult@@@Z PROC ; btSoftSoftCollisionAlgorithm::calculateTimeOfImpact, COMDAT
; _this$ = ecx
; Line 45
	push	ebp
	mov	ebp, esp
	push	ecx
	mov	DWORD PTR [ebp-4], -858993460		; ccccccccH
	mov	DWORD PTR _this$[ebp], ecx
; Line 47
	fld1
; Line 48
	mov	esp, ebp
	pop	ebp
	ret	16					; 00000010H
?calculateTimeOfImpact@btSoftSoftCollisionAlgorithm@@UAEMPAVbtCollisionObject@@0ABUbtDispatcherInfo@@PAVbtManifoldResult@@@Z ENDP ; btSoftSoftCollisionAlgorithm::calculateTimeOfImpact
_TEXT	ENDS
PUBLIC	??2@YAPAXIPAX@Z					; operator new
PUBLIC	?reserve@?$btAlignedObjectArray@PAVbtPersistentManifold@@@@QAEXH@Z ; btAlignedObjectArray<btPersistentManifold *>::reserve
PUBLIC	?allocSize@?$btAlignedObjectArray@PAVbtPersistentManifold@@@@IAEHH@Z ; btAlignedObjectArray<btPersistentManifold *>::allocSize
PUBLIC	?capacity@?$btAlignedObjectArray@PAVbtPersistentManifold@@@@QBEHXZ ; btAlignedObjectArray<btPersistentManifold *>::capacity
PUBLIC	?size@?$btAlignedObjectArray@PAVbtPersistentManifold@@@@QBEHXZ ; btAlignedObjectArray<btPersistentManifold *>::size
; Function compile flags: /Odtp /RTCsu
;	COMDAT ?push_back@?$btAlignedObjectArray@PAVbtPersistentManifold@@@@QAEXABQAVbtPersistentManifold@@@Z
_TEXT	SEGMENT
tv84 = -16						; size = 4
$T20892 = -12						; size = 4
_sz$ = -8						; size = 4
_this$ = -4						; size = 4
__Val$ = 8						; size = 4
?push_back@?$btAlignedObjectArray@PAVbtPersistentManifold@@@@QAEXABQAVbtPersistentManifold@@@Z PROC ; btAlignedObjectArray<btPersistentManifold *>::push_back, COMDAT
; _this$ = ecx
; File d:\專題\自建專案\自建專案\src\linearmath\btalignedobjectarray.h
; Line 277
	push	ebp
	mov	ebp, esp
	sub	esp, 16					; 00000010H
	mov	eax, -858993460				; ccccccccH
	mov	DWORD PTR [ebp-16], eax
	mov	DWORD PTR [ebp-12], eax
	mov	DWORD PTR [ebp-8], eax
	mov	DWORD PTR [ebp-4], eax
	mov	DWORD PTR _this$[ebp], ecx
; Line 278
	mov	ecx, DWORD PTR _this$[ebp]
	call	?size@?$btAlignedObjectArray@PAVbtPersistentManifold@@@@QBEHXZ ; btAlignedObjectArray<btPersistentManifold *>::size
	mov	DWORD PTR _sz$[ebp], eax
; Line 279
	mov	ecx, DWORD PTR _this$[ebp]
	call	?capacity@?$btAlignedObjectArray@PAVbtPersistentManifold@@@@QBEHXZ ; btAlignedObjectArray<btPersistentManifold *>::capacity
	cmp	DWORD PTR _sz$[ebp], eax
	jne	SHORT $LN1@push_back
; Line 281
	mov	ecx, DWORD PTR _this$[ebp]
	call	?size@?$btAlignedObjectArray@PAVbtPersistentManifold@@@@QBEHXZ ; btAlignedObjectArray<btPersistentManifold *>::size
	push	eax
	mov	ecx, DWORD PTR _this$[ebp]
	call	?allocSize@?$btAlignedObjectArray@PAVbtPersistentManifold@@@@IAEHH@Z ; btAlignedObjectArray<btPersistentManifold *>::allocSize
	push	eax
	mov	ecx, DWORD PTR _this$[ebp]
	call	?reserve@?$btAlignedObjectArray@PAVbtPersistentManifold@@@@QAEXH@Z ; btAlignedObjectArray<btPersistentManifold *>::reserve
$LN1@push_back:
; Line 285
	mov	eax, DWORD PTR _this$[ebp]
	mov	ecx, DWORD PTR [eax+4]
	mov	edx, DWORD PTR _this$[ebp]
	mov	eax, DWORD PTR [edx+12]
	lea	ecx, DWORD PTR [eax+ecx*4]
	push	ecx
	push	4
	call	??2@YAPAXIPAX@Z				; operator new
	add	esp, 8
	mov	DWORD PTR $T20892[ebp], eax
	cmp	DWORD PTR $T20892[ebp], 0
	je	SHORT $LN4@push_back
	mov	edx, DWORD PTR $T20892[ebp]
	mov	eax, DWORD PTR __Val$[ebp]
	mov	ecx, DWORD PTR [eax]
	mov	DWORD PTR [edx], ecx
	mov	edx, DWORD PTR $T20892[ebp]
	mov	DWORD PTR tv84[ebp], edx
	jmp	SHORT $LN5@push_back
$LN4@push_back:
	mov	DWORD PTR tv84[ebp], 0
$LN5@push_back:
; Line 290
	mov	eax, DWORD PTR _this$[ebp]
	mov	ecx, DWORD PTR [eax+4]
	add	ecx, 1
	mov	edx, DWORD PTR _this$[ebp]
	mov	DWORD PTR [edx+4], ecx
; Line 291
	add	esp, 16					; 00000010H
	cmp	ebp, esp
	call	__RTC_CheckEsp
	mov	esp, ebp
	pop	ebp
	ret	4
?push_back@?$btAlignedObjectArray@PAVbtPersistentManifold@@@@QAEXABQAVbtPersistentManifold@@@Z ENDP ; btAlignedObjectArray<btPersistentManifold *>::push_back
; Function compile flags: /Odtp /RTCsu
_TEXT	ENDS
;	COMDAT ??2@YAPAXIPAX@Z
_TEXT	SEGMENT
___formal$ = 8						; size = 4
__Where$ = 12						; size = 4
??2@YAPAXIPAX@Z PROC					; operator new, COMDAT
; File c:\program files (x86)\microsoft visual studio 10.0\vc\include\new
; Line 56
	push	ebp
	mov	ebp, esp
; Line 57
	mov	eax, DWORD PTR __Where$[ebp]
; Line 58
	pop	ebp
	ret	0
??2@YAPAXIPAX@Z ENDP					; operator new
; Function compile flags: /Odtp /RTCsu
_TEXT	ENDS
;	COMDAT ?size@?$btAlignedObjectArray@PAVbtPersistentManifold@@@@QBEHXZ
_TEXT	SEGMENT
_this$ = -4						; size = 4
?size@?$btAlignedObjectArray@PAVbtPersistentManifold@@@@QBEHXZ PROC ; btAlignedObjectArray<btPersistentManifold *>::size, COMDAT
; _this$ = ecx
; File d:\專題\自建專案\自建專案\src\linearmath\btalignedobjectarray.h
; Line 150
	push	ebp
	mov	ebp, esp
	push	ecx
	mov	DWORD PTR [ebp-4], -858993460		; ccccccccH
	mov	DWORD PTR _this$[ebp], ecx
; Line 151
	mov	eax, DWORD PTR _this$[ebp]
	mov	eax, DWORD PTR [eax+4]
; Line 152
	mov	esp, ebp
	pop	ebp
	ret	0
?size@?$btAlignedObjectArray@PAVbtPersistentManifold@@@@QBEHXZ ENDP ; btAlignedObjectArray<btPersistentManifold *>::size
; Function compile flags: /Odtp /RTCsu
_TEXT	ENDS
;	COMDAT ?allocSize@?$btAlignedObjectArray@PAVbtPersistentManifold@@@@IAEHH@Z
_TEXT	SEGMENT
tv66 = -8						; size = 4
_this$ = -4						; size = 4
_size$ = 8						; size = 4
?allocSize@?$btAlignedObjectArray@PAVbtPersistentManifold@@@@IAEHH@Z PROC ; btAlignedObjectArray<btPersistentManifold *>::allocSize, COMDAT
; _this$ = ecx
; Line 71
	push	ebp
	mov	ebp, esp
	sub	esp, 8
	mov	DWORD PTR [ebp-8], -858993460		; ccccccccH
	mov	DWORD PTR [ebp-4], -858993460		; ccccccccH
	mov	DWORD PTR _this$[ebp], ecx
; Line 72
	cmp	DWORD PTR _size$[ebp], 0
	je	SHORT $LN3@allocSize
	mov	eax, DWORD PTR _size$[ebp]
	shl	eax, 1
	mov	DWORD PTR tv66[ebp], eax
	jmp	SHORT $LN4@allocSize
$LN3@allocSize:
	mov	DWORD PTR tv66[ebp], 1
$LN4@allocSize:
	mov	eax, DWORD PTR tv66[ebp]
; Line 73
	mov	esp, ebp
	pop	ebp
	ret	4
?allocSize@?$btAlignedObjectArray@PAVbtPersistentManifold@@@@IAEHH@Z ENDP ; btAlignedObjectArray<btPersistentManifold *>::allocSize
; Function compile flags: /Odtp /RTCsu
_TEXT	ENDS
;	COMDAT ?capacity@?$btAlignedObjectArray@PAVbtPersistentManifold@@@@QBEHXZ
_TEXT	SEGMENT
_this$ = -4						; size = 4
?capacity@?$btAlignedObjectArray@PAVbtPersistentManifold@@@@QBEHXZ PROC ; btAlignedObjectArray<btPersistentManifold *>::capacity, COMDAT
; _this$ = ecx
; Line 296
	push	ebp
	mov	ebp, esp
	push	ecx
	mov	DWORD PTR [ebp-4], -858993460		; ccccccccH
	mov	DWORD PTR _this$[ebp], ecx
; Line 297
	mov	eax, DWORD PTR _this$[ebp]
	mov	eax, DWORD PTR [eax+8]
; Line 298
	mov	esp, ebp
	pop	ebp
	ret	0
?capacity@?$btAlignedObjectArray@PAVbtPersistentManifold@@@@QBEHXZ ENDP ; btAlignedObjectArray<btPersistentManifold *>::capacity
_TEXT	ENDS
PUBLIC	?deallocate@?$btAlignedObjectArray@PAVbtPersistentManifold@@@@IAEXXZ ; btAlignedObjectArray<btPersistentManifold *>::deallocate
PUBLIC	?destroy@?$btAlignedObjectArray@PAVbtPersistentManifold@@@@IAEXHH@Z ; btAlignedObjectArray<btPersistentManifold *>::destroy
PUBLIC	?copy@?$btAlignedObjectArray@PAVbtPersistentManifold@@@@IBEXHHPAPAVbtPersistentManifold@@@Z ; btAlignedObjectArray<btPersistentManifold *>::copy
PUBLIC	?allocate@?$btAlignedObjectArray@PAVbtPersistentManifold@@@@IAEPAXH@Z ; btAlignedObjectArray<btPersistentManifold *>::allocate
; Function compile flags: /Odtp /RTCsu
;	COMDAT ?reserve@?$btAlignedObjectArray@PAVbtPersistentManifold@@@@QAEXH@Z
_TEXT	SEGMENT
_s$19971 = -8						; size = 4
_this$ = -4						; size = 4
__Count$ = 8						; size = 4
?reserve@?$btAlignedObjectArray@PAVbtPersistentManifold@@@@QAEXH@Z PROC ; btAlignedObjectArray<btPersistentManifold *>::reserve, COMDAT
; _this$ = ecx
; Line 301
	push	ebp
	mov	ebp, esp
	sub	esp, 8
	mov	DWORD PTR [ebp-8], -858993460		; ccccccccH
	mov	DWORD PTR [ebp-4], -858993460		; ccccccccH
	mov	DWORD PTR _this$[ebp], ecx
; Line 302
	mov	ecx, DWORD PTR _this$[ebp]
	call	?capacity@?$btAlignedObjectArray@PAVbtPersistentManifold@@@@QBEHXZ ; btAlignedObjectArray<btPersistentManifold *>::capacity
	cmp	eax, DWORD PTR __Count$[ebp]
	jge	SHORT $LN2@reserve
; Line 304
	mov	eax, DWORD PTR __Count$[ebp]
	push	eax
	mov	ecx, DWORD PTR _this$[ebp]
	call	?allocate@?$btAlignedObjectArray@PAVbtPersistentManifold@@@@IAEPAXH@Z ; btAlignedObjectArray<btPersistentManifold *>::allocate
	mov	DWORD PTR _s$19971[ebp], eax
; Line 306
	mov	ecx, DWORD PTR _s$19971[ebp]
	push	ecx
	mov	ecx, DWORD PTR _this$[ebp]
	call	?size@?$btAlignedObjectArray@PAVbtPersistentManifold@@@@QBEHXZ ; btAlignedObjectArray<btPersistentManifold *>::size
	push	eax
	push	0
	mov	ecx, DWORD PTR _this$[ebp]
	call	?copy@?$btAlignedObjectArray@PAVbtPersistentManifold@@@@IBEXHHPAPAVbtPersistentManifold@@@Z ; btAlignedObjectArray<btPersistentManifold *>::copy
; Line 308
	mov	ecx, DWORD PTR _this$[ebp]
	call	?size@?$btAlignedObjectArray@PAVbtPersistentManifold@@@@QBEHXZ ; btAlignedObjectArray<btPersistentManifold *>::size
	push	eax
	push	0
	mov	ecx, DWORD PTR _this$[ebp]
	call	?destroy@?$btAlignedObjectArray@PAVbtPersistentManifold@@@@IAEXHH@Z ; btAlignedObjectArray<btPersistentManifold *>::destroy
; Line 310
	mov	ecx, DWORD PTR _this$[ebp]
	call	?deallocate@?$btAlignedObjectArray@PAVbtPersistentManifold@@@@IAEXXZ ; btAlignedObjectArray<btPersistentManifold *>::deallocate
; Line 313
	mov	edx, DWORD PTR _this$[ebp]
	mov	BYTE PTR [edx+16], 1
; Line 315
	mov	eax, DWORD PTR _this$[ebp]
	mov	ecx, DWORD PTR _s$19971[ebp]
	mov	DWORD PTR [eax+12], ecx
; Line 317
	mov	edx, DWORD PTR _this$[ebp]
	mov	eax, DWORD PTR __Count$[ebp]
	mov	DWORD PTR [edx+8], eax
$LN2@reserve:
; Line 320
	add	esp, 8
	cmp	ebp, esp
	call	__RTC_CheckEsp
	mov	esp, ebp
	pop	ebp
	ret	4
?reserve@?$btAlignedObjectArray@PAVbtPersistentManifold@@@@QAEXH@Z ENDP ; btAlignedObjectArray<btPersistentManifold *>::reserve
; Function compile flags: /Odtp /RTCsu
_TEXT	ENDS
;	COMDAT ?copy@?$btAlignedObjectArray@PAVbtPersistentManifold@@@@IBEXHHPAPAVbtPersistentManifold@@@Z
_TEXT	SEGMENT
tv75 = -16						; size = 4
$T20909 = -12						; size = 4
_i$ = -8						; size = 4
_this$ = -4						; size = 4
_start$ = 8						; size = 4
_end$ = 12						; size = 4
_dest$ = 16						; size = 4
?copy@?$btAlignedObjectArray@PAVbtPersistentManifold@@@@IBEXHHPAPAVbtPersistentManifold@@@Z PROC ; btAlignedObjectArray<btPersistentManifold *>::copy, COMDAT
; _this$ = ecx
; Line 75
	push	ebp
	mov	ebp, esp
	sub	esp, 16					; 00000010H
	mov	eax, -858993460				; ccccccccH
	mov	DWORD PTR [ebp-16], eax
	mov	DWORD PTR [ebp-12], eax
	mov	DWORD PTR [ebp-8], eax
	mov	DWORD PTR [ebp-4], eax
	mov	DWORD PTR _this$[ebp], ecx
; Line 77
	mov	eax, DWORD PTR _start$[ebp]
	mov	DWORD PTR _i$[ebp], eax
	jmp	SHORT $LN3@copy
$LN2@copy:
	mov	ecx, DWORD PTR _i$[ebp]
	add	ecx, 1
	mov	DWORD PTR _i$[ebp], ecx
$LN3@copy:
	mov	edx, DWORD PTR _i$[ebp]
	cmp	edx, DWORD PTR _end$[ebp]
	jge	SHORT $LN4@copy
; Line 79
	mov	eax, DWORD PTR _i$[ebp]
	mov	ecx, DWORD PTR _dest$[ebp]
	lea	edx, DWORD PTR [ecx+eax*4]
	push	edx
	push	4
	call	??2@YAPAXIPAX@Z				; operator new
	add	esp, 8
	mov	DWORD PTR $T20909[ebp], eax
	cmp	DWORD PTR $T20909[ebp], 0
	je	SHORT $LN6@copy
	mov	eax, DWORD PTR _this$[ebp]
	mov	ecx, DWORD PTR [eax+12]
	mov	edx, DWORD PTR $T20909[ebp]
	mov	eax, DWORD PTR _i$[ebp]
	mov	ecx, DWORD PTR [ecx+eax*4]
	mov	DWORD PTR [edx], ecx
	mov	edx, DWORD PTR $T20909[ebp]
	mov	DWORD PTR tv75[ebp], edx
	jmp	SHORT $LN7@copy
$LN6@copy:
	mov	DWORD PTR tv75[ebp], 0
$LN7@copy:
	jmp	SHORT $LN2@copy
$LN4@copy:
; Line 83
	add	esp, 16					; 00000010H
	cmp	ebp, esp
	call	__RTC_CheckEsp
	mov	esp, ebp
	pop	ebp
	ret	12					; 0000000cH
?copy@?$btAlignedObjectArray@PAVbtPersistentManifold@@@@IBEXHHPAPAVbtPersistentManifold@@@Z ENDP ; btAlignedObjectArray<btPersistentManifold *>::copy
; Function compile flags: /Odtp /RTCsu
_TEXT	ENDS
;	COMDAT ?destroy@?$btAlignedObjectArray@PAVbtPersistentManifold@@@@IAEXHH@Z
_TEXT	SEGMENT
_i$ = -8						; size = 4
_this$ = -4						; size = 4
_first$ = 8						; size = 4
_last$ = 12						; size = 4
?destroy@?$btAlignedObjectArray@PAVbtPersistentManifold@@@@IAEXHH@Z PROC ; btAlignedObjectArray<btPersistentManifold *>::destroy, COMDAT
; _this$ = ecx
; Line 94
	push	ebp
	mov	ebp, esp
	sub	esp, 8
	mov	DWORD PTR [ebp-8], -858993460		; ccccccccH
	mov	DWORD PTR [ebp-4], -858993460		; ccccccccH
	mov	DWORD PTR _this$[ebp], ecx
; Line 96
	mov	eax, DWORD PTR _first$[ebp]
	mov	DWORD PTR _i$[ebp], eax
	jmp	SHORT $LN3@destroy
$LN2@destroy:
	mov	ecx, DWORD PTR _i$[ebp]
	add	ecx, 1
	mov	DWORD PTR _i$[ebp], ecx
$LN3@destroy:
	mov	edx, DWORD PTR _i$[ebp]
	cmp	edx, DWORD PTR _last$[ebp]
	jge	SHORT $LN4@destroy
; Line 99
	jmp	SHORT $LN2@destroy
$LN4@destroy:
; Line 100
	mov	esp, ebp
	pop	ebp
	ret	8
?destroy@?$btAlignedObjectArray@PAVbtPersistentManifold@@@@IAEXHH@Z ENDP ; btAlignedObjectArray<btPersistentManifold *>::destroy
_TEXT	ENDS
PUBLIC	?allocate@?$btAlignedAllocator@PAVbtPersistentManifold@@$0BA@@@QAEPAPAVbtPersistentManifold@@HPAPBQAV2@@Z ; btAlignedAllocator<btPersistentManifold *,16>::allocate
; Function compile flags: /Odtp /RTCsu
;	COMDAT ?allocate@?$btAlignedObjectArray@PAVbtPersistentManifold@@@@IAEPAXH@Z
_TEXT	SEGMENT
_this$ = -4						; size = 4
_size$ = 8						; size = 4
?allocate@?$btAlignedObjectArray@PAVbtPersistentManifold@@@@IAEPAXH@Z PROC ; btAlignedObjectArray<btPersistentManifold *>::allocate, COMDAT
; _this$ = ecx
; Line 103
	push	ebp
	mov	ebp, esp
	push	ecx
	mov	DWORD PTR [ebp-4], -858993460		; ccccccccH
	mov	DWORD PTR _this$[ebp], ecx
; Line 104
	cmp	DWORD PTR _size$[ebp], 0
	je	SHORT $LN1@allocate
; Line 105
	push	0
	mov	eax, DWORD PTR _size$[ebp]
	push	eax
	mov	ecx, DWORD PTR _this$[ebp]
	call	?allocate@?$btAlignedAllocator@PAVbtPersistentManifold@@$0BA@@@QAEPAPAVbtPersistentManifold@@HPAPBQAV2@@Z ; btAlignedAllocator<btPersistentManifold *,16>::allocate
	jmp	SHORT $LN2@allocate
$LN1@allocate:
; Line 106
	xor	eax, eax
$LN2@allocate:
; Line 107
	add	esp, 4
	cmp	ebp, esp
	call	__RTC_CheckEsp
	mov	esp, ebp
	pop	ebp
	ret	4
?allocate@?$btAlignedObjectArray@PAVbtPersistentManifold@@@@IAEPAXH@Z ENDP ; btAlignedObjectArray<btPersistentManifold *>::allocate
_TEXT	ENDS
PUBLIC	?deallocate@?$btAlignedAllocator@PAVbtPersistentManifold@@$0BA@@@QAEXPAPAVbtPersistentManifold@@@Z ; btAlignedAllocator<btPersistentManifold *,16>::deallocate
; Function compile flags: /Odtp /RTCsu
;	COMDAT ?deallocate@?$btAlignedObjectArray@PAVbtPersistentManifold@@@@IAEXXZ
_TEXT	SEGMENT
_this$ = -4						; size = 4
?deallocate@?$btAlignedObjectArray@PAVbtPersistentManifold@@@@IAEXXZ PROC ; btAlignedObjectArray<btPersistentManifold *>::deallocate, COMDAT
; _this$ = ecx
; Line 110
	push	ebp
	mov	ebp, esp
	push	ecx
	mov	DWORD PTR [ebp-4], -858993460		; ccccccccH
	mov	DWORD PTR _this$[ebp], ecx
; Line 111
	mov	eax, DWORD PTR _this$[ebp]
	cmp	DWORD PTR [eax+12], 0
	je	SHORT $LN3@deallocate
; Line 113
	mov	ecx, DWORD PTR _this$[ebp]
	movzx	edx, BYTE PTR [ecx+16]
	test	edx, edx
	je	SHORT $LN1@deallocate
; Line 115
	mov	eax, DWORD PTR _this$[ebp]
	mov	ecx, DWORD PTR [eax+12]
	push	ecx
	mov	ecx, DWORD PTR _this$[ebp]
	call	?deallocate@?$btAlignedAllocator@PAVbtPersistentManifold@@$0BA@@@QAEXPAPAVbtPersistentManifold@@@Z ; btAlignedAllocator<btPersistentManifold *,16>::deallocate
$LN1@deallocate:
; Line 117
	mov	edx, DWORD PTR _this$[ebp]
	mov	DWORD PTR [edx+12], 0
$LN3@deallocate:
; Line 119
	add	esp, 4
	cmp	ebp, esp
	call	__RTC_CheckEsp
	mov	esp, ebp
	pop	ebp
	ret	0
?deallocate@?$btAlignedObjectArray@PAVbtPersistentManifold@@@@IAEXXZ ENDP ; btAlignedObjectArray<btPersistentManifold *>::deallocate
_TEXT	ENDS
EXTRN	?btAlignedAllocInternal@@YAPAXIH@Z:PROC		; btAlignedAllocInternal
; Function compile flags: /Odtp /RTCsu
;	COMDAT ?allocate@?$btAlignedAllocator@PAVbtPersistentManifold@@$0BA@@@QAEPAPAVbtPersistentManifold@@HPAPBQAV2@@Z
_TEXT	SEGMENT
_this$ = -4						; size = 4
_n$ = 8							; size = 4
_hint$ = 12						; size = 4
?allocate@?$btAlignedAllocator@PAVbtPersistentManifold@@$0BA@@@QAEPAPAVbtPersistentManifold@@HPAPBQAV2@@Z PROC ; btAlignedAllocator<btPersistentManifold *,16>::allocate, COMDAT
; _this$ = ecx
; File d:\專題\自建專案\自建專案\src\linearmath\btalignedallocator.h
; Line 84
	push	ebp
	mov	ebp, esp
	push	ecx
	mov	DWORD PTR [ebp-4], -858993460		; ccccccccH
	mov	DWORD PTR _this$[ebp], ecx
; Line 86
	push	16					; 00000010H
	mov	eax, DWORD PTR _n$[ebp]
	shl	eax, 2
	push	eax
	call	?btAlignedAllocInternal@@YAPAXIH@Z	; btAlignedAllocInternal
	add	esp, 8
; Line 87
	add	esp, 4
	cmp	ebp, esp
	call	__RTC_CheckEsp
	mov	esp, ebp
	pop	ebp
	ret	8
?allocate@?$btAlignedAllocator@PAVbtPersistentManifold@@$0BA@@@QAEPAPAVbtPersistentManifold@@HPAPBQAV2@@Z ENDP ; btAlignedAllocator<btPersistentManifold *,16>::allocate
_TEXT	ENDS
EXTRN	?btAlignedFreeInternal@@YAXPAX@Z:PROC		; btAlignedFreeInternal
; Function compile flags: /Odtp /RTCsu
;	COMDAT ?deallocate@?$btAlignedAllocator@PAVbtPersistentManifold@@$0BA@@@QAEXPAPAVbtPersistentManifold@@@Z
_TEXT	SEGMENT
_this$ = -4						; size = 4
_ptr$ = 8						; size = 4
?deallocate@?$btAlignedAllocator@PAVbtPersistentManifold@@$0BA@@@QAEXPAPAVbtPersistentManifold@@@Z PROC ; btAlignedAllocator<btPersistentManifold *,16>::deallocate, COMDAT
; _this$ = ecx
; Line 89
	push	ebp
	mov	ebp, esp
	push	ecx
	mov	DWORD PTR [ebp-4], -858993460		; ccccccccH
	mov	DWORD PTR _this$[ebp], ecx
; Line 90
	mov	eax, DWORD PTR _ptr$[ebp]
	push	eax
	call	?btAlignedFreeInternal@@YAXPAX@Z	; btAlignedFreeInternal
	add	esp, 4
; Line 91
	add	esp, 4
	cmp	ebp, esp
	call	__RTC_CheckEsp
	mov	esp, ebp
	pop	ebp
	ret	4
?deallocate@?$btAlignedAllocator@PAVbtPersistentManifold@@$0BA@@@QAEXPAPAVbtPersistentManifold@@@Z ENDP ; btAlignedAllocator<btPersistentManifold *,16>::deallocate
_TEXT	ENDS
END
