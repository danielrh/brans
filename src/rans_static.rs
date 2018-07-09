#![feature(core_intrinsics)]
type c_void = u8;
type c_char = u8;
type c_uchar = u8;
type c_schar = u8;
type c_ulong = u64;
type c_long = i64;
type c_ushort = u16;
type c_short = i16;
type c_uint = u32;
type c_int = i32;
type c_double = f64;
type c_float = f32;
extern crate libc;
extern crate core;
extern "C" {
    #[no_mangle]
    fn malloc(_: c_ulong) -> *mut c_void;
    #[no_mangle]
    fn memmove(_: *mut c_void, _: *const c_void, _: c_ulong)
     -> *mut c_void;
    #[no_mangle]
    fn memset(_: *mut c_void, _: c_int, _: c_ulong)
     -> *mut c_void;
}
pub type uint32_t = c_uint;
pub type __off_t = c_long;
pub type Rans64State = uint64_t;

#[derive ( Copy , Clone )]
#[repr ( C )]
pub struct ari_decoder {
    pub R: [c_uchar; 4096],
}
#[derive ( Copy , Clone , Debug)]
#[repr ( C )]
pub struct RansEncSymbol {
    pub rcp_freq: uint64_t,
    pub freq: uint32_t,
    pub bias: uint32_t,
    pub cmpl_freq: uint32_t,
    pub rcp_shift: uint32_t,
}
pub type Rans64DecSymbol = RansDecSymbol;
#[derive ( Copy , Clone, Debug)]
#[repr ( C )]
pub struct RansDecSymbol {
    pub start: u16,
    pub freq: u16,
}

#[derive ( Copy , Clone, Debug)]
pub struct RansDecSymbolStartFreq {
    pub start: u16,
    pub freq: u16,
    pub sym: u8,    
}

pub type RansState = Rans64State;
pub type Rans64EncSymbol = RansEncSymbol;
pub type uint64_t = c_ulong;
pub type uint8_t = c_uchar;
#[derive ( Copy , Clone )]
#[repr ( C )]
pub struct blocks {
    pub blk: *mut c_uchar,
    pub sz: uint32_t,
}
fn Rans64MulHi(mut a: uint64_t, mut b: uint64_t)
 -> uint64_t {
    return ((a as u128).wrapping_mul(b as u128) >> 64i32) as uint64_t;
}
#[no_mangle]
fn Rans64EncInit(mut r: &mut Rans64State) -> () {
    *r = (1u64 << 31i32) as Rans64State;
}
#[no_mangle]
fn Rans64EncPut(mut r: &mut Rans64State,
                                  mut pptr: &mut [uint32_t],
                                  mut ptr_off: &mut usize, // pre-decrement
                                  mut start: uint32_t, mut freq: uint32_t,
                                  mut scale_bits: uint32_t) -> () {
    let mut x: uint64_t = *r;
    let mut x_max: uint64_t =
        (1u64 << 31i32 >> scale_bits <<
             32i32).wrapping_mul(freq as u64) as uint64_t;
    if x >= x_max {
        //eprintln!("E:{:x}", x as uint32_t);
        pptr[*ptr_off] = x as uint32_t;
        *ptr_off = ptr_off.wrapping_sub(1);
        x >>= 32i32
    }
    *r =
        (x.wrapping_div(freq as c_ulong) <<
             scale_bits).wrapping_add(x.wrapping_rem(freq as
                                                         c_ulong)).wrapping_add(start
                                                                                          as
                                                                                          c_ulong);
}
#[no_mangle]
fn Rans64EncFlush(mut r: &mut Rans64State,
                                    mut pptr: &mut [uint32_t],
                                    ptr_off: &mut usize) {
    let mut x: uint64_t = *r;
    pptr[ptr_off.wrapping_sub(1)] = x as uint32_t;
    pptr[*ptr_off] = (x >> 32) as uint32_t;
    *ptr_off -= 2;
    //eprintln!("V:{:x}", x);
}
#[no_mangle]
fn Rans64DecInit(mut r: &mut Rans64State,
                 pptr: &[uint32_t],
                 ptr_off: &mut usize) -> () {
    let mut x: uint64_t = 0;
    x = pptr[*ptr_off] as u64;
    x |= (pptr[*ptr_off + 1] as u64) << 32;
    *ptr_off += 2;
    *r = x;
}
#[no_mangle]
fn Rans64DecGet(mut r: &mut Rans64State,
                mut scale_bits: uint32_t) -> uint32_t {
    return (*r &
                (1u32 << scale_bits).wrapping_sub(1i32 as c_uint) as
                    c_ulong) as uint32_t;
}
#[no_mangle]
fn Rans64DecAdvance(mut r: &mut Rans64State,
                    pptr: &[uint32_t],
                    ptr_off: &mut usize,
                    mut start: u16, mut freq: u16,
                    mut scale_bits: uint32_t) -> () {
    let mut mask: uint64_t =
        (1u64 << scale_bits).wrapping_sub(1i32 as u64) as
            uint64_t;
    let mut x: uint64_t = *r;
    x =
        (freq as
             c_ulong).wrapping_mul(x >>
                                             scale_bits).wrapping_add(x &
                                                                          mask).wrapping_sub(start
                                                                                                 as
                                                                                                 c_ulong);
    if (x as u64) < 1u64 << 31i32 {
        x = x << 32i32 | pptr[*ptr_off] as c_ulong;
        *ptr_off = ptr_off.wrapping_add(1);
    }
    *r = x;
}



///////////////u8/////////
#[no_mangle]
fn u8Rans64DecInit(mut r: &mut Rans64State,
                 pptr: &[u8],
                 ptr_off: &mut u32) -> () {
    let mut x: uint64_t = 0;
    x = (pptr[*ptr_off as usize + 7] as u64) << 56;
    x |= (pptr[*ptr_off as usize + 6] as u64) << 48;
    x |= (pptr[*ptr_off as usize + 5] as u64) << 40;
    x |= (pptr[*ptr_off as usize + 4] as u64) << 32;
    x |= (pptr[*ptr_off as usize + 3] as u64) << 24;
    x |= (pptr[*ptr_off as usize + 2] as u64) << 16;
    x |= (pptr[*ptr_off as usize + 1] as u64) << 8;
    x |= (pptr[*ptr_off as usize] as u64);
    //eprintln!("I{:x}", x);
    *ptr_off = ptr_off.wrapping_add(8);
    *r = x;
}

#[no_mangle]
fn u8Rans64DecAdvance(mut r: &mut Rans64State,
                      pptr: &[u8],
                      ptr_off: &mut u32,
                      mut start: u16, mut freq: u16,
                      mut scale_bits: uint32_t) -> () {
    let mut mask: uint64_t =
        (1u64 << scale_bits).wrapping_sub(1i32 as u64) as
            uint64_t;
    let mut x: uint64_t = *r;
    x =
        (freq as
             c_ulong).wrapping_mul(x >>
                                             scale_bits).wrapping_add(x &
                                                                          mask).wrapping_sub(start
                                                                                                 as
                                                                                                 c_ulong);
    if (x as u64) < 1u64 << 31i32 {
        x = x << 32i32 | (
            (pptr[*ptr_off as usize + 3] as u64) << 24) | (
            (pptr[*ptr_off as usize + 2] as u64) << 16) | (
            (pptr[*ptr_off as usize + 1] as u64) << 8) | (
            pptr[*ptr_off as usize] as u64);
        //eprintln!("D{:x}", x &0xffff_ffff);
        *ptr_off = ptr_off.wrapping_add(4);
    }
    *r = x;
}

#[no_mangle]
fn u8Rans64DecAdvanceSymbol(mut r: &mut Rans64State,
                                            pptr: &[u8],
                                            ptr_off: &mut u32,
                                            mut sym: &Rans64DecSymbol,
                                            mut scale_bits: uint32_t) {
    u8Rans64DecAdvance(r, pptr, ptr_off, (*sym).start, (*sym).freq, scale_bits);
}

#[no_mangle]
fn u8Rans64DecRenorm(mut r: &mut Rans64State,
                   pptr: &[u8], ptr_offset: &mut u32) {
    let mut x: uint64_t = *r;
    if (x as u64) < 1u64 << 31i32 {
        x = x << 32i32 | (
            (pptr[*ptr_offset as usize + 3] as u64) << 24) | (
            (pptr[*ptr_offset as usize + 2] as u64) << 16) | (
            (pptr[*ptr_offset as usize + 1] as u64) << 8) | (
            pptr[*ptr_offset as usize] as u64);
        //eprintln!("D{:x}", x &0xffff_ffff);
        *ptr_offset = ptr_offset.wrapping_add(4);
    }
    *r = x;
}
#[no_mangle]
fn u8Rans64DecForceRenorm(mut r: &mut Rans64State,
                   pptr: &[u8], ptr_offset: &mut u32) {
    let mut x: uint64_t = *r;
    x = x << 32i32 | (
        (pptr[*ptr_offset as usize + 3] as u64) << 24) | (
        (pptr[*ptr_offset as usize + 2] as u64) << 16) | (
        (pptr[*ptr_offset as usize + 1] as u64) << 8) | (
        pptr[*ptr_offset as usize] as u64);
    //eprintln!("D{:x}", x &0xffff_ffff);
    *ptr_offset = ptr_offset.wrapping_add(4);
    *r = x;
}
///////////////u8///////////


















#[no_mangle]
fn Rans64EncSymbolInit(mut s: &mut Rans64EncSymbol,
                                         mut start: uint32_t,
                                         mut freq: uint32_t,
                                         mut scale_bits: uint32_t) -> () {
    (*s).freq = freq;
    (*s).cmpl_freq =
        ((1i32 << scale_bits) as c_uint).wrapping_sub(freq);
    if freq < 2i32 as c_uint {
        (*s).rcp_freq = !0u64 as uint64_t;
        (*s).rcp_shift = 0i32 as uint32_t;
        (*s).bias =
            start.wrapping_add((1i32 << scale_bits) as
                                   c_uint).wrapping_sub(1i32 as
                                                                  c_uint)
    } else {
        let mut shift: uint32_t = 0i32 as uint32_t;
        let mut x0: uint64_t = 0;
        let mut x1: uint64_t = 0;
        let mut t0: uint64_t = 0;
        let mut t1: uint64_t = 0;
        while freq > 1u32 << shift { shift = shift.wrapping_add(1) }
        x0 = freq.wrapping_sub(1i32 as c_uint) as uint64_t;
        x1 = (1u64 << shift.wrapping_add(31i32 as c_uint)) as uint64_t;
        t1 = x1.wrapping_div(freq as c_ulong);
        x0 =
            (x0 as
                 c_ulong).wrapping_add(x1.wrapping_rem(freq as
                                                                 c_ulong)
                                                 << 32i32) as uint64_t as
                uint64_t;
        t0 = x0.wrapping_div(freq as c_ulong);
        (*s).rcp_freq = t0.wrapping_add(t1 << 32i32);
        (*s).rcp_shift = shift.wrapping_sub(1i32 as c_uint);
        (*s).bias = start
    };
}
#[no_mangle]
fn Rans64DecSymbolInit(mut s: &mut Rans64DecSymbol,
                                         mut start: u16,
                                         mut freq: u16) -> () {
    (*s).start = start;
    (*s).freq = freq;
}

#[no_mangle]
fn Rans64EncPutSymbol(mut r: &mut Rans64State,
                                        pptr: &mut [uint32_t],
                                        mut ptr_offset: &mut usize,
                                        mut sym: &Rans64EncSymbol,
                                        mut scale_bits: uint32_t) -> () {
    let mut x: uint64_t = *r;
    let mut x_max: uint64_t =
        (1u64 << 31i32 >> scale_bits <<
             32i32).wrapping_mul((*sym).freq as u64) as
            uint64_t;
    if x >= x_max {
        pptr[*ptr_offset] = x as u32;
        //eprintln!("E:{:x}", x as uint32_t);
        *ptr_offset = ptr_offset.wrapping_sub(1);
        x >>= 32i32
    }
    let mut q: uint64_t = Rans64MulHi(x, (*sym).rcp_freq) >> (*sym).rcp_shift;
    *r =
        x.wrapping_add((*sym).bias as
                           c_ulong).wrapping_add(q.wrapping_mul((*sym).cmpl_freq
                                                                          as
                                                                          c_ulong));
}
#[no_mangle]
fn Rans64DecAdvanceSymbol(mut r: &mut Rans64State,
                                            pptr: &[uint32_t],
                                            ptr_off: &mut usize,
                                            mut sym: &Rans64DecSymbol,
                                            mut scale_bits: uint32_t) -> () {
    Rans64DecAdvance(r, pptr, ptr_off, (*sym).start, (*sym).freq, scale_bits);
}
fn Rans64DecAdvanceStep(mut r: &mut Rans64State,
                                          mut start: u16,
                                          mut freq: u16,
                                          mut scale_bits: uint32_t) -> () {
    let mut mask: uint64_t =
        (1u32 << scale_bits).wrapping_sub(1i32 as c_uint) as uint64_t;
    let mut x: uint64_t = *r;
    *r =
        (freq as
             c_ulong).wrapping_mul(x >>
                                             scale_bits).wrapping_add(x &
                                                                          mask).wrapping_sub(start
                                                                                                 as
                                                                                                 c_ulong);
}
#[no_mangle]
fn Rans64DecAdvanceSymbolStep(mut r: &mut Rans64State,
                                                mut sym:
                                                    &Rans64DecSymbol,
                                                mut scale_bits: uint32_t)
 -> () {
    Rans64DecAdvanceStep(r, (*sym).start, (*sym).freq, scale_bits);
}
#[no_mangle]
fn Rans64DecRenorm(mut r: &mut Rans64State,
                   pptr: &[uint32_t], ptr_offset: &mut usize) -> () {
    let mut x: uint64_t = *r;
    if (x as u64) < 1u64 << 31i32 {
        x = x << 32i32 | pptr[*ptr_offset] as c_ulong;
        *ptr_offset = ptr_offset.wrapping_add(1);
    }
    *r = x;
}
#[no_mangle]
unsafe extern "C" fn hist8(mut in_0: *mut c_uchar,
                           mut in_size: c_uint,
                           mut F0: *mut c_int) -> () {
    let mut F1: [c_int; 264] = [0;264];
    let mut F2: [c_int; 264] = [0;264];
    let mut F3: [c_int; 264] = [0;264];
    let mut F4: [c_int; 264] = [0;264];
    let mut F5: [c_int; 264] = [0;264];
    let mut F6: [c_int; 264] = [0;264];
    let mut F7: [c_int; 264] = [0;264];
    let mut i: c_int = 0;
    let mut i8: c_int =
        (in_size & !8i32 as c_uint) as c_int;
    i = 0i32;
    while i < i8 {
        let ref mut fresh10 =
            *F0.offset(*in_0.offset((i + 0i32) as isize) as isize);
        *fresh10 += 1;
        F1[*in_0.offset((i + 1i32) as isize) as usize] += 1;
        F2[*in_0.offset((i + 2i32) as isize) as usize] += 1;
        F3[*in_0.offset((i + 3i32) as isize) as usize] += 1;
        F4[*in_0.offset((i + 4i32) as isize) as usize] += 1;
        F5[*in_0.offset((i + 5i32) as isize) as usize] += 1;
        F6[*in_0.offset((i + 6i32) as isize) as usize] += 1;
        F7[*in_0.offset((i + 7i32) as isize) as usize] += 1;
        i += 8i32
    }
    while (i as c_uint) < in_size {
        let fresh11 = i;
        i = i + 1;
        let ref mut fresh12 =
            *F0.offset(*in_0.offset(fresh11 as isize) as isize);
        *fresh12 += 1
    }
    i = 0i32;
    while i < 256i32 {
        *F0.offset(i as isize) +=
            F1[i as usize] + F2[i as usize] + F3[i as usize] + F4[i as usize]
                + F5[i as usize] + F6[i as usize] + F7[i as usize];
        i += 1
    };
}
#[no_mangle]
pub unsafe extern "C" fn rans_compress_O0(mut in_0: *mut c_uchar,
                                          mut in_size: c_uint,
                                          mut out_size: *mut c_uint)
 -> *mut c_uchar {
    let mut current_block: u64;
    let mut sz: c_int =
        (1.05f64 * in_size as c_double +
             (257i32 * 257i32 * 3i32) as c_double +
             8i32 as c_double) as c_int & !3i32;
    let mut out_buf: *mut c_uchar =
        malloc((sz + 4i32) as c_ulong) as *mut c_uchar;
    let mut cp: *mut c_uchar = 0 as *mut c_uchar;
    let mut out_end: *mut c_uchar = 0 as *mut c_uchar;
    let mut syms: [RansEncSymbol; 256] =
        [RansEncSymbol{rcp_freq: 0,
                       freq: 0,
                       bias: 0,
                       cmpl_freq: 0,
                       rcp_shift: 0,}; 256];
    let mut rans0: RansState = 0;
    let mut rans1: RansState = 0;
    let mut rans2: RansState = 0;
    let mut rans3: RansState = 0;
    let mut ptr: *mut uint8_t = 0 as *mut uint8_t;
    let mut F: [c_int; 264] = [0;264];
    let mut i: c_int = 0;
    let mut j: c_int = 0;
    let mut tab_size: c_int = 0;
    let mut rle: c_int = 0;
    let mut x: c_int = 0;
    let mut fsum: c_int = 0i32;
    let mut m: c_int = 0i32;
    let mut M: c_int = 0i32;
    let mut tr: uint64_t = 0;
    if out_buf.is_null() {
        return 0 as *mut c_uchar
    } else {
        out_end = out_buf.offset(sz as isize);
        ptr = out_end;
        hist8(in_0, in_size, F.as_mut_ptr());
        tr =
            (((1i32 << 12i32) as uint64_t) <<
                 31i32).wrapping_div(in_size as
                                         c_ulong).wrapping_add(((1i32 <<
                                                                           30i32)
                                                                          as
                                                                          c_uint).wrapping_div(in_size)
                                                                         as
                                                                         c_ulong);
        j = 0i32;
        M = j;
        m = M;
        while j < 256i32 {
            if !(0 == F[j as usize]) {
                if m < F[j as usize] { m = F[j as usize]; M = j }
                F[j as usize] =
                    ((F[j as usize] as c_ulong).wrapping_mul(tr) >>
                         31i32) as c_int;
                if F[j as usize] == 0i32 { F[j as usize] = 1i32 }
                fsum += F[j as usize]
            }
            j += 1
        }
        fsum += 1;
        if fsum < 1i32 << 12i32 {
            F[M as usize] += (1i32 << 12i32) - fsum
        } else { F[M as usize] -= fsum - (1i32 << 12i32) }
        cp = out_buf.offset(4isize);
        j = 0i32;
        rle = j;
        x = rle;
        while j < 256i32 {
            if 0 != F[j as usize] {
                if 0 != rle {
                    rle -= 1
                } else {
                    let fresh13 = cp;
                    cp = cp.offset(1);
                    *fresh13 = j as c_uchar;
                    if 0 == rle && 0 != j && 0 != F[(j - 1i32) as usize] {
                        rle = j + 1i32;
                        while rle < 256i32 && 0 != F[rle as usize] {
                            rle += 1
                        }
                        rle -= j + 1i32;
                        let fresh14 = cp;
                        cp = cp.offset(1);
                        *fresh14 = rle as c_uchar
                    }
                }
                if F[j as usize] < 128i32 {
                    let fresh15 = cp;
                    cp = cp.offset(1);
                    *fresh15 = F[j as usize] as c_uchar
                } else {
                    let fresh16 = cp;
                    cp = cp.offset(1);
                    *fresh16 =
                        (128i32 | F[j as usize] >> 8i32) as c_uchar;
                    let fresh17 = cp;
                    cp = cp.offset(1);
                    *fresh17 = (F[j as usize] & 255i32) as c_uchar
                }
                Rans64EncSymbolInit(&mut syms[j as usize]
                                      , x as uint32_t,
                                    F[j as usize] as uint32_t,
                                    12i32 as uint32_t);
                x += F[j as usize]
            }
            j += 1
        }
        let fresh18 = cp;
        cp = cp.offset(1);
        *fresh18 = 0i32 as c_uchar;
        tab_size =
            out_buf.offset_to(cp).expect("bad offset_to") as c_long as
                c_int;
        Rans64EncInit(&mut rans0);
        Rans64EncInit(&mut rans1);
        Rans64EncInit(&mut rans2);
        Rans64EncInit(&mut rans3);
        i = (in_size & 3i32 as c_uint) as c_int;
        match i {
            3 => {
                let mut new_offset=  0usize;
                let pptr = core::slice::from_raw_parts_mut((ptr as *mut uint8_t as *mut uint32_t).offset(-1), 1);
                Rans64EncPutSymbol(&mut rans2,
                                   pptr, &mut new_offset,
                                   &mut syms[*in_0.offset(in_size.wrapping_sub((i
                                                                                    -
                                                                                    2i32)
                                                                                   as
                                                                                   c_uint)
                                                              as isize) as
                                                 usize],
                                   12i32 as uint32_t);
                ptr = (ptr as *mut uint8_t as *mut uint32_t).offset(new_offset as isize) as *mut c_uchar;
                current_block = 14453171571987214953;
            }
            2 => { current_block = 14453171571987214953; }
            1 => { current_block = 15409834050583242653; }
            0 | _ => { current_block = 224731115979188411; }
        }
        match current_block {
            14453171571987214953 => {
                let mut new_offset=  0usize;
                let pptr = core::slice::from_raw_parts_mut((ptr as *mut uint8_t as *mut uint32_t).offset(-1), 1);
                Rans64EncPutSymbol(&mut rans1,
                                   pptr, &mut new_offset,
                                   &mut syms[*in_0.offset(in_size.wrapping_sub((i
                                                                                    -
                                                                                    1i32)
                                                                                   as
                                                                                   c_uint)
                                                              as isize) as
                                                 usize],
                                   12i32 as uint32_t);
                ptr = (ptr as *mut uint8_t as *mut uint32_t).offset(new_offset as isize) as *mut c_uchar;
                current_block = 15409834050583242653;
            }
            _ => { }
        }
        match current_block {
            15409834050583242653 => {
                let mut new_offset=  0usize;
                let pptr = core::slice::from_raw_parts_mut((ptr as *mut uint8_t as *mut uint32_t).offset(-1), 1);
                Rans64EncPutSymbol(&mut rans0,
                                   pptr, &mut new_offset,
                                   &mut syms[*in_0.offset(in_size.wrapping_sub((i
                                                                                    -
                                                                                    0i32)
                                                                                   as
                                                                                   c_uint)
                                                              as isize) as
                                                 usize],
                                   12i32 as uint32_t);
                ptr = (ptr as *mut uint8_t as *mut uint32_t).offset(new_offset as isize) as *mut c_uchar;
            }
            _ => { }
        }
        i = (in_size & !3i32 as c_uint) as c_int;
        while i > 0i32 {
            let mut new_offset=  3usize;
            let pptr = core::slice::from_raw_parts_mut((ptr as *mut uint8_t as *mut uint32_t).offset(-4), 4);
            {
                let mut s3: &mut RansEncSymbol =
                    &mut syms[*in_0.offset((i - 1i32) as isize) as usize];
                Rans64EncPutSymbol(&mut rans3,
                                   pptr, &mut new_offset,
                                   s3, 12i32 as uint32_t);
            }
            {
                let mut s2: &mut RansEncSymbol =
                    &mut syms[*in_0.offset((i - 2i32) as isize) as usize];
                Rans64EncPutSymbol(&mut rans2,
                                   pptr, &mut new_offset,
                                   s2, 12i32 as uint32_t);
            }
            {
                let mut s1: &mut RansEncSymbol =
                    &mut syms[*in_0.offset((i - 3i32) as isize) as usize];
                Rans64EncPutSymbol(&mut rans1,
                                   pptr, &mut new_offset,
                                   s1, 12i32 as uint32_t);
            }
            let mut s0: &mut RansEncSymbol =
                &mut syms[*in_0.offset((i - 4i32) as isize) as usize];
            Rans64EncPutSymbol(&mut rans0,
                               pptr, &mut new_offset,
                               s0, 12i32 as uint32_t);
            ptr = (ptr as *mut uint8_t as *mut uint32_t).offset(new_offset as isize - 3) as *mut c_uchar;
            i -= 4i32
        }
        let mut new_offset=  7usize;
        let pptr = core::slice::from_raw_parts_mut((ptr as *mut uint8_t as *mut uint32_t).offset(-8), 4);
        Rans64EncFlush(&mut rans3,
                       pptr, &mut new_offset);
        Rans64EncFlush(&mut rans2,
                       pptr, &mut new_offset);
        Rans64EncFlush(&mut rans1,
                       pptr, &mut new_offset);
        Rans64EncFlush(&mut rans0,
                       pptr, &mut new_offset);
        ptr = (ptr as *mut uint8_t as *mut uint32_t).offset(new_offset as isize - 3) as *mut c_uchar;
        
        *out_size =
            (ptr.offset_to(out_end).expect("bad offset_to") as c_long +
                 tab_size as c_long) as c_uint;
        cp = out_buf;
        let fresh19 = cp;
        cp = cp.offset(1);
        *fresh19 =
            (in_size >> 0i32 & 255i32 as c_uint) as c_uchar;
        let fresh20 = cp;
        cp = cp.offset(1);
        *fresh20 =
            (in_size >> 8i32 & 255i32 as c_uint) as c_uchar;
        let fresh21 = cp;
        cp = cp.offset(1);
        *fresh21 =
            (in_size >> 16i32 & 255i32 as c_uint) as c_uchar;
        let fresh22 = cp;
        cp = cp.offset(1);
        *fresh22 =
            (in_size >> 24i32 & 255i32 as c_uint) as c_uchar;
        memmove(out_buf.offset(tab_size as isize) as *mut c_void,
                ptr as *const c_void,
                ptr.offset_to(out_end).expect("bad offset_to") as c_long
                    as c_ulong);
        return out_buf
    };
}
/*
#[no_mangle]
pub unsafe extern "C" fn rans_uncompress_O0(mut in_0: *mut c_uchar,
                                            mut in_size: c_uint,
                                            mut out_size: *mut c_uint)
 -> *mut c_uchar {
    let mut c_0: c_uchar = 0;
    let mut cp: *mut c_uchar = in_0.offset(4isize);
    let mut ht_in_offset = 4usize;
    let mut i: c_int = 0;
    let mut j: c_int = 0;
    let mut x: c_int = 0;
    let mut out_sz: c_int = 0;
    let mut rle: c_int = 0;
    let mut out_buf: *mut c_char = 0 as *mut c_char;
    let mut D: ari_decoder = ari_decoder{R: [0; 4096],};
    let mut syms: [RansDecSymbol; 256] =
        [RansDecSymbol{start: 0, freq: 0,}; 256];
    out_sz =
        (*in_0.offset(0isize) as c_int) << 0i32 |
            (*in_0.offset(1isize) as c_int) << 8i32 |
            (*in_0.offset(2isize) as c_int) << 16i32 |
            (*in_0.offset(3isize) as c_int) << 24i32;
    out_buf = malloc(out_sz as c_ulong) as *mut c_char;
    if out_buf.is_null() {
        return 0 as *mut c_uchar
    } else {
        x = 0i32;
        rle = x;
        let fresh23 = cp;
        cp = cp.offset(1);
        ht_in_offset += 1;
        j = *fresh23 as c_int;
        loop  {
            let mut F: c_int = 0;
            let mut C: c_int = 0;
            let fresh24 = cp;
            cp = cp.offset(1);
            ht_in_offset += 1;
            F = *fresh24 as c_int;
            if F >= 128i32 {
                F &= !128i32;
                let fresh25 = cp;
                cp = cp.offset(1);
                ht_in_offset += 1;
                F = (F & 127i32) << 8i32 | *fresh25 as c_int
            }
            C = x;
            Rans64DecSymbolInit(&mut syms[j as usize],
                                C as uint32_t, F as uint32_t);
            memset(&mut D.R[x as usize] as *mut c_uchar as
                       *mut c_void, j, F as c_ulong);
            x += F;
            if 0 == rle && j + 1i32 == *cp as c_int {
                let fresh26 = cp;
                cp = cp.offset(1);
                ht_in_offset += 1;
                j = *fresh26 as c_int;
                let fresh27 = cp;
                cp = cp.offset(1);
                ht_in_offset += 1;
                rle = *fresh27 as c_int
            } else if 0 != rle {
                rle -= 1;
                j += 1
            } else {
                let fresh28 = cp;
                cp = cp.offset(1);
                ht_in_offset += 1;
                j = *fresh28 as c_int
            }
            if !(0 != j) { break ; }
        }
        let mut rans0: RansState = 0;
        let mut rans1: RansState = 0;
        let mut rans2: RansState = 0;
        let mut rans3: RansState = 0;
        let mut ptr: *mut uint8_t = cp;
        let pptr = core::slice::from_raw_parts_mut(ptr, in_size - ht_in_offset);
        let ptr_offset = 0usize;

        Rans64DecInit(&mut rans0,
                      pptr, &mut ptr_offset);
        Rans64DecInit(&mut rans1,
                      pptr, &mut ptr_offset);
        Rans64DecInit(&mut rans2,
                      pptr, &mut ptr_offset);
        Rans64DecInit(&mut rans3,
                      pptr, &mut ptr_offset);
        let mut out_end: c_int = out_sz & !3i32;
        let mut R: [RansState; 4] = [0; 4];
        R[0usize] = rans0;
        R[1usize] = rans1;
        R[2usize] = rans2;
        R[3usize] = rans3;
        let mut mask: uint32_t =
            (1u32 << 12i32).wrapping_sub(1i32 as c_uint);
        i = 0i32;
        while i < out_end {
            let mut m: [uint32_t; 4] = [0; 4];
            let mut c: [uint8_t; 4] = [0; 4];
            m[0usize] = (R[0usize] & mask as c_ulong) as uint32_t;
            c[0usize] = D.R[m[0usize] as usize];
            *out_buf.offset((i + 0i32) as isize) = c[0usize] as c_char;
            R[0usize] =
                (syms[c[0usize] as usize].freq as
                     c_ulong).wrapping_mul(R[0usize] >> 12i32);
            R[0usize] =
                (R[0usize] as
                     c_ulong).wrapping_add(m[0usize].wrapping_sub(syms[c[0usize]
                                                                                 as
                                                                                 usize].start)
                                                     as c_ulong) as
                    RansState as RansState;
            m[1usize] = (R[1usize] & mask as c_ulong) as uint32_t;
            c[1usize] = D.R[m[1usize] as usize];
            *out_buf.offset((i + 1i32) as isize) = c[1usize] as c_char;
            R[1usize] =
                (syms[c[1usize] as usize].freq as
                     c_ulong).wrapping_mul(R[1usize] >> 12i32);
            R[1usize] =
                (R[1usize] as
                     c_ulong).wrapping_add(m[1usize].wrapping_sub(syms[c[1usize]
                                                                                 as
                                                                                 usize].start)
                                                     as c_ulong) as
                    RansState as RansState;
            m[2usize] = (R[2usize] & mask as c_ulong) as uint32_t;
            c[2usize] = D.R[m[2usize] as usize];
            *out_buf.offset((i + 2i32) as isize) = c[2usize] as c_char;
            R[2usize] =
                (syms[c[2usize] as usize].freq as
                     c_ulong).wrapping_mul(R[2usize] >> 12i32);
            R[2usize] =
                (R[2usize] as
                     c_ulong).wrapping_add(m[2usize].wrapping_sub(syms[c[2usize]
                                                                                 as
                                                                                 usize].start)
                                                     as c_ulong) as
                    RansState as RansState;
            m[3usize] = (R[3usize] & mask as c_ulong) as uint32_t;
            c[3usize] = D.R[m[3usize] as usize];
            *out_buf.offset((i + 3i32) as isize) = c[3usize] as c_char;
            R[3usize] =
                (syms[c[3usize] as usize].freq as
                     c_ulong).wrapping_mul(R[3usize] >> 12i32);
            R[3usize] =
                (R[3usize] as
                     c_ulong).wrapping_add(m[3usize].wrapping_sub(syms[c[3usize]
                                                                                 as
                                                                                 usize].start)
                                                     as c_ulong) as
                    RansState as RansState;
            Rans64DecRenorm(&mut R[0usize],
                            pptr, &mut ptr_offset);
            Rans64DecRenorm(&mut R[1usize],
                            pptr, &mut ptr_offset);
            Rans64DecRenorm(&mut R[2usize],
                            pptr, &mut ptr_offset);
            Rans64DecRenorm(&mut R[3usize],
                            pptr, &mut ptr_offset);
            i += 4i32
        }
        rans0 = R[0usize];
        rans1 = R[1usize];
        rans2 = R[2usize];
        rans3 = R[3usize];
        match out_sz & 3i32 {
            1 => {
                c_0 =
                    D.R[Rans64DecGet(&mut rans0,
                                     12i32 as uint32_t) as usize];
                Rans64DecAdvanceSymbol(&mut rans0,
                                       pptr, &mut ptr_offset,
                                       &mut syms[c_0 as usize],
                                       12i32 as uint32_t);
                *out_buf.offset(out_end as isize) = c_0 as c_char
            }
            2 => {
                c_0 =
                    D.R[Rans64DecGet(&mut rans0,
                                     12i32 as uint32_t) as usize];
                Rans64DecAdvanceSymbol(&mut rans0,
                                       pptr, &mut ptr_offset,
                                       &mut syms[c_0 as usize],
                                       12i32 as uint32_t);
                *out_buf.offset(out_end as isize) = c_0 as c_char;
                c_0 =
                    D.R[Rans64DecGet(&mut rans1,
                                     12i32 as uint32_t) as usize];
                Rans64DecAdvanceSymbol(&mut rans1,
                                       pptr, &mut ptr_offset,
                                       &mut syms[c_0 as usize],
                                       12i32 as uint32_t);
                *out_buf.offset((out_end + 1i32) as isize) =
                    c_0 as c_char
            }
            3 => {
                c_0 =
                    D.R[Rans64DecGet(&mut rans0,
                                     12i32 as uint32_t) as usize];
                Rans64DecAdvanceSymbol(&mut rans0,
                                       pptr, &mut ptr_offset,
                                       &mut syms[c_0 as usize],
                                       12i32 as uint32_t);
                *out_buf.offset(out_end as isize) = c_0 as c_char;
                c_0 =
                    D.R[Rans64DecGet(&mut rans1,
                                     12i32 as uint32_t) as usize];
                Rans64DecAdvanceSymbol(&mut rans1,
                                       pptr, &mut ptr_offset,
                                       &mut syms[c_0 as usize],
                                       12i32 as uint32_t);
                *out_buf.offset((out_end + 1i32) as isize) =
                    c_0 as c_char;
                c_0 =
                    D.R[Rans64DecGet(&mut rans2,
                                     12i32 as uint32_t) as usize];
                Rans64DecAdvanceSymbol(&mut rans2,
                                       pptr, &mut ptr_offset,
                                       &mut syms[c_0 as usize],
                                       12i32 as uint32_t);
                *out_buf.offset((out_end + 2i32) as isize) =
                    c_0 as c_char
            }
            0 | _ => { }
        }
        *out_size = out_sz as c_uint;
        return out_buf as *mut c_uchar
    };
}
*/
#[no_mangle]
pub unsafe extern "C" fn rans_uncompress_O0b(mut in_0: *mut c_uchar,
                                             mut in_size: c_uint,
                                             mut out_size: *mut c_uint)
 -> *mut c_uchar {
    let mut c_0: c_uchar = 0;
    let mut cp: *mut c_uchar = in_0.offset(4isize);
    let mut ht_in_offset = 4;
    let mut i: c_int = 0;
    let mut j: c_int = 0;
    let mut x: c_int = 0;
    let mut out_sz: c_int = 0;
    let mut rle: c_int = 0;
    let mut out_buf: *mut c_char = 0 as *mut c_char;
    let mut D: ari_decoder = ari_decoder{R: [0; 4096],};
    let mut syms: [RansDecSymbol; 256] =
        [RansDecSymbol{start: 0, freq: 0,}; 256];
    out_sz =
        (*in_0.offset(0isize) as c_int) << 0i32 |
            (*in_0.offset(1isize) as c_int) << 8i32 |
            (*in_0.offset(2isize) as c_int) << 16i32 |
            (*in_0.offset(3isize) as c_int) << 24i32;
    out_buf = malloc(out_sz as c_ulong) as *mut c_char;
    if out_buf.is_null() {
        return 0 as *mut c_uchar
    } else {
        x = 0i32;
        rle = x;
        let fresh29 = cp;
        cp = cp.offset(1);
        ht_in_offset += 1;
        j = *fresh29 as c_int;
        loop  {
            let mut F: c_int = 0;
            let mut C: c_int = 0;
            let fresh30 = cp;
            cp = cp.offset(1);
            ht_in_offset += 1;
            F = *fresh30 as c_int;
            if F >= 128i32 {
                F &= !128i32;
                let fresh31 = cp;
                cp = cp.offset(1);
                ht_in_offset += 1;
                F = (F & 127i32) << 8i32 | *fresh31 as c_int
            }
            C = x;
            Rans64DecSymbolInit(&mut syms[j as usize],
                                C as u16, F as u16);
            memset(&mut D.R[x as usize] as *mut c_uchar as
                       *mut c_void, j, F as c_ulong);
            x += F;
            if 0 == rle && j + 1i32 == *cp as c_int {
                let fresh33 = cp;
                cp = cp.offset(1);
                ht_in_offset += 1;
                j = *fresh33 as c_int;
                let fresh34 = cp;
                cp = cp.offset(1);
                ht_in_offset += 1;
                rle = *fresh34 as c_int
            } else if 0 != rle {
                rle -= 1;
                j += 1
            } else {
                let fresh35 = cp;
                cp = cp.offset(1);
                ht_in_offset += 1;
                j = *fresh35 as c_int
            }
            if !(0 != j) { break ; }
        }
        let mut rans0: RansState = 0;
        let mut rans1: RansState = 0;
        let mut rans2: RansState = 0;
        let mut rans3: RansState = 0;
        let mut ptr: *mut uint8_t = cp;
        let pptr = core::slice::from_raw_parts_mut(ptr as *mut u32, (in_size as usize - ht_in_offset) >> 2);
        let mut ptr_offset = 0usize;

        Rans64DecInit(&mut rans0,
                      pptr, &mut ptr_offset);
        Rans64DecInit(&mut rans1,
                      pptr, &mut ptr_offset);
        
        Rans64DecInit(&mut rans2,
                      pptr, &mut ptr_offset);
        Rans64DecInit(&mut rans3,
                      pptr, &mut ptr_offset);
        let mut out_end: c_int = out_sz & !3i32;
        let mut R: [RansState; 4] = [0; 4];
        R[0usize] = rans0;
        R[1usize] = rans1;
        R[2usize] = rans2;
        R[3usize] = rans3;
        let mut mask: uint32_t =
            (1u32 << 12i32).wrapping_sub(1i32 as c_uint);
        i = 0i32;
        while i < out_end {
            let mut m: [uint32_t; 4] =
                [(R[0usize] & mask as c_ulong) as uint32_t,
                 (R[1usize] & mask as c_ulong) as uint32_t,
                 (R[2usize] & mask as c_ulong) as uint32_t,
                 (R[3usize] & mask as c_ulong) as uint32_t];
            let mut c: [uint8_t; 4] =
                [D.R[m[0usize] as usize], D.R[m[1usize] as usize],
                 D.R[m[2usize] as usize], D.R[m[3usize] as usize]];
            *out_buf.offset((i + 0i32) as isize) = c[0usize] as c_char;
            *out_buf.offset((i + 1i32) as isize) = c[1usize] as c_char;
            *out_buf.offset((i + 2i32) as isize) = c[2usize] as c_char;
            *out_buf.offset((i + 3i32) as isize) = c[3usize] as c_char;
            R[0usize] =
                (syms[c[0usize] as usize].freq as
                     c_ulong).wrapping_mul(R[0usize] >> 12i32);
            R[0usize] =
                (R[0usize] as
                     c_ulong).wrapping_add(m[0usize].wrapping_sub(u32::from(syms[c[0usize]
                                                                                 as
                                                                                 usize].start))
                                                     as c_ulong) as
                    RansState as RansState;
            R[1usize] =
                (syms[c[1usize] as usize].freq as
                     c_ulong).wrapping_mul(R[1usize] >> 12i32);
            R[1usize] =
                (R[1usize] as
                     c_ulong).wrapping_add(m[1usize].wrapping_sub(u32::from(syms[c[1usize]
                                                                                 as
                                                                                 usize].start))
                                                     as c_ulong) as
                    RansState as RansState;
            R[2usize] =
                (syms[c[2usize] as usize].freq as
                     c_ulong).wrapping_mul(R[2usize] >> 12i32);
            R[2usize] =
                (R[2usize] as
                     c_ulong).wrapping_add(m[2usize].wrapping_sub(u32::from(syms[c[2usize]
                                                                                 as
                                                                                 usize].start))
                                                     as c_ulong) as
                    RansState as RansState;
            R[3usize] =
                (syms[c[3usize] as usize].freq as
                     c_ulong).wrapping_mul(R[3usize] >> 12i32);
            R[3usize] =
                (R[3usize] as
                     c_ulong).wrapping_add(m[3usize].wrapping_sub(u32::from(syms[c[3usize]
                                                                                 as
                                                                                 usize].start))
                                                     as c_ulong) as
                    RansState as RansState;
            Rans64DecRenorm(&mut R[0usize],
                            pptr, &mut ptr_offset);
            Rans64DecRenorm(&mut R[1usize],
                            pptr, &mut ptr_offset);
            Rans64DecRenorm(&mut R[2usize],
                            pptr, &mut ptr_offset);
            Rans64DecRenorm(&mut R[3usize],
                            pptr, &mut ptr_offset);

            i += 4i32
        }
        rans0 = R[0usize];
        rans1 = R[1usize];
        rans2 = R[2usize];
        rans3 = R[3usize];
        match out_sz & 3i32 {
            1 => {
                c_0 =
                    D.R[Rans64DecGet(&mut rans0,
                                     12i32 as uint32_t) as usize];
                Rans64DecAdvanceSymbol(&mut rans0,
                                       pptr, &mut ptr_offset,
                                       &mut syms[c_0 as usize],
                                       12i32 as uint32_t);
                *out_buf.offset(out_end as isize) = c_0 as c_char
            }
            2 => {
                c_0 =
                    D.R[Rans64DecGet(&mut rans0,
                                     12i32 as uint32_t) as usize];
                Rans64DecAdvanceSymbol(&mut rans0,
                                       pptr, &mut ptr_offset,
                                       &mut syms[c_0 as usize],
                                       12i32 as uint32_t);
                *out_buf.offset(out_end as isize) = c_0 as c_char;
                c_0 =
                    D.R[Rans64DecGet(&mut rans1,
                                     12i32 as uint32_t) as usize];
                Rans64DecAdvanceSymbol(&mut rans1,
                                       pptr, &mut ptr_offset,
                                       &mut syms[c_0 as usize],
                                       12i32 as uint32_t);
                *out_buf.offset((out_end + 1i32) as isize) =
                    c_0 as c_char
            }
            3 => {
                c_0 =
                    D.R[Rans64DecGet(&mut rans0,
                                     12i32 as uint32_t) as usize];
                Rans64DecAdvanceSymbol(&mut rans0,
                                       pptr, &mut ptr_offset,
                                       &mut syms[c_0 as usize],
                                       12i32 as uint32_t);
                *out_buf.offset(out_end as isize) = c_0 as c_char;
                c_0 =
                    D.R[Rans64DecGet(&mut rans1,
                                     12i32 as uint32_t) as usize];
                Rans64DecAdvanceSymbol(&mut rans1,
                                       pptr, &mut ptr_offset,
                                       &mut syms[c_0 as usize],
                                       12i32 as uint32_t);
                *out_buf.offset((out_end + 1i32) as isize) =
                    c_0 as c_char;
                c_0 =
                    D.R[Rans64DecGet(&mut rans2,
                                     12i32 as uint32_t) as usize];
                Rans64DecAdvanceSymbol(&mut rans2,
                                       pptr, &mut ptr_offset,
                                       &mut syms[c_0 as usize],
                                       12i32 as uint32_t);
                *out_buf.offset((out_end + 2i32) as isize) =
                    c_0 as c_char
            }
            0 | _ => { }
        }
        *out_size = out_sz as c_uint;
        return out_buf as *mut c_uchar
    };
}
unsafe extern "C" fn hist1_4(mut in_0: *mut c_uchar,
                             mut in_size: c_uint,
                             mut F0: *mut [c_int; 256],
                             mut T0: *mut c_int) -> () {
    let mut T1: [c_int; 264] = [0;264];
    let mut T2: [c_int; 264] = [0;264];
    let mut T3: [c_int; 264] = [0;264];
    let mut idiv4: c_uint = in_size.wrapping_div(4i32 as c_uint);
    let mut i: c_int = 0;
    let mut c0: c_uchar = 0;
    let mut c1: c_uchar = 0;
    let mut c2: c_uchar = 0;
    let mut c3: c_uchar = 0;
    let mut in0: *mut c_uchar = in_0.offset(0isize);
    let mut in1: *mut c_uchar = in_0.offset(idiv4 as isize);
    let mut in2: *mut c_uchar =
        in_0.offset(idiv4.wrapping_mul(2i32 as c_uint) as isize);
    let mut in3: *mut c_uchar =
        in_0.offset(idiv4.wrapping_mul(3i32 as c_uint) as isize);
    let mut last_0: c_uchar = 0i32 as c_uchar;
    let mut last_1: c_uchar = *in1.offset(-1i32 as isize);
    let mut last_2: c_uchar = *in2.offset(-1i32 as isize);
    let mut last_3: c_uchar = *in3.offset(-1i32 as isize);
    let mut in0_end: *mut c_uchar = in1;
    while in0 < in0_end {
        let fresh36 = in0;
        in0 = in0.offset(1);
        c0 = *fresh36;
        let ref mut fresh37 = (*F0.offset(last_0 as isize))[c0 as usize];
        *fresh37 += 1;
        let ref mut fresh38 = *T0.offset(last_0 as isize);
        *fresh38 += 1;
        last_0 = c0;
        let fresh39 = in1;
        in1 = in1.offset(1);
        c1 = *fresh39;
        let ref mut fresh40 = (*F0.offset(last_1 as isize))[c1 as usize];
        *fresh40 += 1;
        T1[last_1 as usize] += 1;
        last_1 = c1;
        let fresh41 = in2;
        in2 = in2.offset(1);
        c2 = *fresh41;
        let ref mut fresh42 = (*F0.offset(last_2 as isize))[c2 as usize];
        *fresh42 += 1;
        T2[last_2 as usize] += 1;
        last_2 = c2;
        let fresh43 = in3;
        in3 = in3.offset(1);
        c3 = *fresh43;
        let ref mut fresh44 = (*F0.offset(last_3 as isize))[c3 as usize];
        *fresh44 += 1;
        T3[last_3 as usize] += 1;
        last_3 = c3
    }
    while in3 < in_0.offset(in_size as isize) {
        let fresh45 = in3;
        in3 = in3.offset(1);
        c3 = *fresh45;
        let ref mut fresh46 = (*F0.offset(last_3 as isize))[c3 as usize];
        *fresh46 += 1;
        T3[last_3 as usize] += 1;
        last_3 = c3
    }
    i = 0i32;
    while i < 256i32 {
        *T0.offset(i as isize) +=
            T1[i as usize] + T2[i as usize] + T3[i as usize];
        i += 1
    };
}
#[no_mangle]
pub unsafe extern "C" fn rans_compress_O1(mut in_0: *mut c_uchar,
                                          mut in_size: c_uint,
                                          mut out_size: *mut c_uint)
 -> *mut c_uchar {
    let mut out_buf: *mut c_uchar =
        malloc((1.05f64 * in_size as c_double +
                    (257i32 * 257i32 * 3i32) as c_double +
                    4i32 as c_double) as c_ulong) as
            *mut c_uchar;
    let mut cp: *mut c_uchar = out_buf;
    let mut out_end: *mut c_uchar = 0 as *mut c_uchar;
    let mut tab_size: c_uint = 0;
    let mut rle_i: c_uint = 0;
    let mut rle_j: c_uint = 0;
    let mut out_end_byte_offset = (1.05f64 * in_size as c_double) as usize + (257i32 * 257i32 * 3i32) as usize + 4;
    let mut syms: [[RansEncSymbol; 256]; 256] =
        [[RansEncSymbol{rcp_freq: 0,
                        freq: 0,
                        bias: 0,
                        cmpl_freq: 0,
                        rcp_shift: 0,}; 256]; 256];
    if out_buf.is_null() {
        return 0 as *mut c_uchar
    } else {
        out_end =
            out_buf.offset((1.05f64 * in_size as c_double) as
                               c_int as
                               isize).offset((257i32 * 257i32 * 3i32) as
                                             isize).offset(4isize);
        assert_eq!(out_end, out_buf.offset(out_end_byte_offset as isize));
        cp = out_buf.offset(4isize);
        let mut F: [[c_int; 256]; 256] = [[0;256];256];
        let mut T: [c_int; 264] = [0;264];
        let mut i: c_int = 0;
        let mut j: c_int = 0;
        hist1_4(in_0, in_size, F.as_mut_ptr(), T.as_mut_ptr());
        F[0usize][*in_0.offset((1i32 as
                                    c_uint).wrapping_mul(in_size >>
                                                                   2i32) as
                                   isize) as usize] += 1;
        F[0usize][*in_0.offset((2i32 as
                                    c_uint).wrapping_mul(in_size >>
                                                                   2i32) as
                                   isize) as usize] += 1;
        F[0usize][*in_0.offset((3i32 as
                                    c_uint).wrapping_mul(in_size >>
                                                                   2i32) as
                                   isize) as usize] += 1;
        T[0usize] += 3i32;
        i = 0i32;
        rle_i = i as c_uint;
        while i < 256i32 {
            let mut t2: c_int = 0;
            let mut m: c_int = 0;
            let mut M: c_int = 0;
            let mut x: c_uint = 0;
            if !(T[i as usize] == 0i32) {
                let mut p: c_double =
                    (1i32 << 12i32) as c_double /
                        T[i as usize] as c_double;
                j = 0i32;
                M = j;
                m = M;
                t2 = m;
                while j < 256i32 {
                    if !(0 == F[i as usize][j as usize]) {
                        if m < F[i as usize][j as usize] {
                            m = F[i as usize][j as usize];
                            M = j
                        }
                        F[i as usize][j as usize] =
                            (F[i as usize][j as usize] as c_double * p)
                                as c_int;
                        if F[i as usize][j as usize] == 0i32 {
                            F[i as usize][j as usize] = 1i32
                        }
                        t2 += F[i as usize][j as usize]
                    }
                    j += 1
                }
                t2 += 1;
                if t2 < 1i32 << 12i32 {
                    F[i as usize][M as usize] += (1i32 << 12i32) - t2
                } else { F[i as usize][M as usize] -= t2 - (1i32 << 12i32) }
                if 0 != rle_i {
                    rle_i = rle_i.wrapping_sub(1)
                } else {
                    let fresh47 = cp;
                    cp = cp.offset(1);
                    *fresh47 = i as c_uchar;
                    if 0 != i && 0 != T[(i - 1i32) as usize] {
                        rle_i = (i + 1i32) as c_uint;
                        while rle_i < 256i32 as c_uint &&
                                  0 != T[rle_i as usize] {
                            rle_i = rle_i.wrapping_add(1)
                        }
                        rle_i =
                            rle_i.wrapping_sub((i + 1i32) as c_uint);
                        let fresh48 = cp;
                        cp = cp.offset(1);
                        *fresh48 = rle_i as c_uchar
                    }
                }
                let mut F_i_: *mut c_int = F[i as usize].as_mut_ptr();
                x = 0i32 as c_uint;
                rle_j = 0i32 as c_uint;
                j = 0i32;
                while j < 256i32 {
                    if 0 != *F_i_.offset(j as isize) {
                        if 0 != rle_j {
                            rle_j = rle_j.wrapping_sub(1)
                        } else {
                            let fresh49 = cp;
                            cp = cp.offset(1);
                            *fresh49 = j as c_uchar;
                            if 0 == rle_j && 0 != j &&
                                   0 != *F_i_.offset((j - 1i32) as isize) {
                                rle_j = (j + 1i32) as c_uint;
                                while rle_j < 256i32 as c_uint &&
                                          0 != *F_i_.offset(rle_j as isize) {
                                    rle_j = rle_j.wrapping_add(1)
                                }
                                rle_j =
                                    rle_j.wrapping_sub((j + 1i32) as
                                                           c_uint);
                                let fresh50 = cp;
                                cp = cp.offset(1);
                                *fresh50 = rle_j as c_uchar
                            }
                        }
                        if *F_i_.offset(j as isize) < 128i32 {
                            let fresh51 = cp;
                            cp = cp.offset(1);
                            *fresh51 =
                                *F_i_.offset(j as isize) as c_uchar
                        } else {
                            let fresh52 = cp;
                            cp = cp.offset(1);
                            *fresh52 =
                                (128i32 | *F_i_.offset(j as isize) >> 8i32) as
                                    c_uchar;
                            let fresh53 = cp;
                            cp = cp.offset(1);
                            *fresh53 =
                                (*F_i_.offset(j as isize) & 255i32) as
                                    c_uchar
                        }
                        Rans64EncSymbolInit(&mut syms[i as usize][j as usize]
                                            , x,
                                            *F_i_.offset(j as isize) as
                                                uint32_t, 12i32 as uint32_t);
                        x =
                            x.wrapping_add(*F_i_.offset(j as isize) as
                                               c_uint)
                    }
                    j += 1
                }
                let fresh54 = cp;
                cp = cp.offset(1);
                *fresh54 = 0i32 as c_uchar
            }
            i += 1
        }
        let fresh55 = cp;
        cp = cp.offset(1);
        *fresh55 = 0i32 as c_uchar;
        tab_size =
            out_buf.offset_to(cp).expect("bad offset_to") as c_long as
                c_uint;
        let mut rans0: RansState = 0;
        let mut rans1: RansState = 0;
        let mut rans2: RansState = 0;
        let mut rans3: RansState = 0;
        Rans64EncInit(&mut rans0);
        Rans64EncInit(&mut rans1);
        Rans64EncInit(&mut rans2);
        Rans64EncInit(&mut rans3);
        let mut ptr: *mut uint8_t = out_end;
        let mut isz4: c_int = (in_size >> 2i32) as c_int;
        let mut i0: c_int = 1i32 * isz4 - 2i32;
        let mut i1: c_int = 2i32 * isz4 - 2i32;
        let mut i2: c_int = 3i32 * isz4 - 2i32;
        let mut i3: c_int = 4i32 * isz4 - 2i32;
        let mut l0: c_uchar = *in_0.offset((i0 + 1i32) as isize);
        let mut l1: c_uchar = *in_0.offset((i1 + 1i32) as isize);
        let mut l2: c_uchar = *in_0.offset((i2 + 1i32) as isize);
        let mut l3: c_uchar = *in_0.offset((i3 + 1i32) as isize);
        let pptr = core::slice::from_raw_parts_mut(out_buf.offset(out_end_byte_offset as isize &3) as *mut uint32_t, out_end_byte_offset>>2);
        let mut ptr_offset = (out_end_byte_offset >> 2) - 1;
        l3 =
            *in_0.offset(in_size.wrapping_sub(1i32 as c_uint) as isize);
        i3 = in_size.wrapping_sub(2i32 as c_uint) as c_int;
        while i3 > 4i32 * isz4 - 2i32 {
            let mut c3: c_uchar = *in_0.offset(i3 as isize);
            Rans64EncPutSymbol(&mut rans3,
                               pptr, &mut ptr_offset,
                               &mut syms[c3 as usize][l3 as usize]
                                  , 12i32 as uint32_t);
            l3 = c3;
            i3 -= 1
        }
        while i0 >= 0i32 {
            let mut c0: c_uchar = 0;
            let mut c1: c_uchar = 0;
            let mut c2: c_uchar = 0;
            let mut c3_0: c_uchar = 0;
            {
                c3_0 = *in_0.offset(i3 as isize);
                let mut s3 =
                    &mut syms[c3_0 as usize][l3 as usize];
                Rans64EncPutSymbol(&mut rans3,
                               pptr, &mut ptr_offset,
                               s3, 12i32 as uint32_t);
            }
            {
            c2 = *in_0.offset(i2 as isize);
                let mut s2 =
                    &mut syms[c2 as usize][l2 as usize];
                Rans64EncPutSymbol(&mut rans2,
                               pptr, &mut ptr_offset,
                               s2, 12i32 as uint32_t);
            }
            {
                c1 = *in_0.offset(i1 as isize);
                let mut s1 =
                    &mut syms[c1 as usize][l1 as usize];
                Rans64EncPutSymbol(&mut rans1,
                               pptr, &mut ptr_offset,
                               s1, 12i32 as uint32_t);
            }
            c0 = *in_0.offset(i0 as isize);
            let mut s0 =
                &mut syms[c0 as usize][l0 as usize];
            Rans64EncPutSymbol(&mut rans0,
                               pptr, &mut ptr_offset,
                               s0, 12i32 as uint32_t);
            l0 = c0;
            l1 = c1;
            l2 = c2;
            l3 = c3_0;
            i0 -= 1;
            i1 -= 1;
            i2 -= 1;
            i3 -= 1
        }
        Rans64EncPutSymbol(&mut rans3,
                           pptr, &mut ptr_offset,
                           &mut syms[0usize][l3 as usize], 12i32 as uint32_t);
        Rans64EncPutSymbol(&mut rans2,
                           pptr, &mut ptr_offset,
                           &mut syms[0usize][l2 as usize], 12i32 as uint32_t);
        Rans64EncPutSymbol(&mut rans1,
                           pptr, &mut ptr_offset,
                           &mut syms[0usize][l1 as usize], 12i32 as uint32_t);
        Rans64EncPutSymbol(&mut rans0,
                           pptr, &mut ptr_offset,
                           &mut syms[0usize][l0 as usize], 12i32 as uint32_t);
        Rans64EncFlush(&mut rans3,
                       pptr, &mut ptr_offset);
        Rans64EncFlush(&mut rans2,
                       pptr, &mut ptr_offset);
        Rans64EncFlush(&mut rans1,
                       pptr, &mut ptr_offset);
        Rans64EncFlush(&mut rans0,
                       pptr, &mut ptr_offset);
        let out_byte_size = (pptr.len() - (ptr_offset.wrapping_add(1))) << 2;
        *out_size = out_byte_size as u32 + tab_size;
        cp = out_buf;
        let fresh56 = cp;
        cp = cp.offset(1);
        *fresh56 =
            (in_size >> 0i32 & 255i32 as c_uint) as c_uchar;
        let fresh57 = cp;
        cp = cp.offset(1);
        *fresh57 =
            (in_size >> 8i32 & 255i32 as c_uint) as c_uchar;
        let fresh58 = cp;
        cp = cp.offset(1);
        *fresh58 =
            (in_size >> 16i32 & 255i32 as c_uint) as c_uchar;
        let fresh59 = cp;
        cp = cp.offset(1);
        *fresh59 =
            (in_size >> 24i32 & 255i32 as c_uint) as c_uchar;
        memmove(out_buf.offset(tab_size as isize) as *mut c_void,
                pptr.as_ptr().offset(ptr_offset as isize + 1) as *const c_void,
                *out_size as c_ulong);
        return out_buf
    };
}
const THREE:usize = 3;

#[no_mangle]
pub unsafe extern "C" fn rans_compress_O1d(mut in_0: *mut c_uchar,
                                          mut in_size: c_uint,
                                          mut out_size: *mut c_uint)
 -> *mut c_uchar {
    let mut out_buf: *mut c_uchar =
        malloc((1.05f64 * in_size as c_double +
                    (257i32 * 257i32 * 3i32) as c_double +
                    4i32 as c_double) as c_ulong) as
            *mut c_uchar;
    let mut cp: *mut c_uchar = out_buf;
    let mut out_end: *mut c_uchar = 0 as *mut c_uchar;
    let mut tab_size: c_uint = 0;
    let mut rle_i: c_uint = 0;
    let mut rle_j: c_uint = 0;
    let mut out_end_byte_offset = (1.05f64 * in_size as c_double) as usize + (257i32 * 257i32 * 3i32) as usize + 4;
    let mut syms: [[RansEncSymbol; 256]; 256] =
        [[RansEncSymbol{rcp_freq: 0,
                        freq: 0,
                        bias: 0,
                        cmpl_freq: 0,
                        rcp_shift: 0,}; 256]; 256];
    if out_buf.is_null() {
        return 0 as *mut c_uchar
    } else {
        out_end =
            out_buf.offset((1.05f64 * in_size as c_double) as
                               c_int as
                               isize).offset((257i32 * 257i32 * 3i32) as
                                             isize).offset(4isize);
        assert_eq!(out_end, out_buf.offset(out_end_byte_offset as isize));
        cp = out_buf.offset(4isize);
        let mut F: [[c_int; 256]; 256] = [[0;256];256];
        let mut T: [c_int; 264] = [0;264];
        let mut i: c_int = 0;
        let mut j: c_int = 0;
        hist1_4(in_0, in_size, F.as_mut_ptr(), T.as_mut_ptr());
        F[0usize][*in_0.offset((1i32 as
                                    c_uint).wrapping_mul(in_size >>
                                                                   2i32) as
                                   isize) as usize] += 1;
        F[0usize][*in_0.offset((2i32 as
                                    c_uint).wrapping_mul(in_size >>
                                                                   2i32) as
                                   isize) as usize] += 1;
        F[0usize][*in_0.offset((3i32 as
                                    c_uint).wrapping_mul(in_size >>
                                                                   2i32) as
                                   isize) as usize] += 1;
        T[0usize] += 3i32;
        i = 0i32;
        rle_i = i as c_uint;
        while i < 256i32 {
            let mut t2: c_int = 0;
            let mut m: c_int = 0;
            let mut M: c_int = 0;
            let mut x: c_uint = 0;
            if !(T[i as usize] == 0i32) {
                let mut p: c_double =
                    (1i32 << 12i32) as c_double /
                        T[i as usize] as c_double;
                j = 0i32;
                M = j;
                m = M;
                t2 = m;
                while j < 256i32 {
                    if !(0 == F[i as usize][j as usize]) {
                        if m < F[i as usize][j as usize] {
                            m = F[i as usize][j as usize];
                            M = j
                        }
                        F[i as usize][j as usize] =
                            (F[i as usize][j as usize] as c_double * p)
                                as c_int;
                        if F[i as usize][j as usize] == 0i32 {
                            F[i as usize][j as usize] = 1i32
                        }
                        t2 += F[i as usize][j as usize]
                    }
                    j += 1
                }
                t2 += 1;
                if t2 < 1i32 << 12i32 {
                    F[i as usize][M as usize] += (1i32 << 12i32) - t2
                } else { F[i as usize][M as usize] -= t2 - (1i32 << 12i32) }
                if 0 != rle_i {
                    rle_i = rle_i.wrapping_sub(1)
                } else {
                    let fresh47 = cp;
                    cp = cp.offset(1);
                    *fresh47 = i as c_uchar;
                    if 0 != i && 0 != T[(i - 1i32) as usize] {
                        rle_i = (i + 1i32) as c_uint;
                        while rle_i < 256i32 as c_uint &&
                                  0 != T[rle_i as usize] {
                            rle_i = rle_i.wrapping_add(1)
                        }
                        rle_i =
                            rle_i.wrapping_sub((i + 1i32) as c_uint);
                        let fresh48 = cp;
                        cp = cp.offset(1);
                        *fresh48 = rle_i as c_uchar
                    }
                }
                let mut F_i_: *mut c_int = F[i as usize].as_mut_ptr();
                x = 0i32 as c_uint;
                rle_j = 0i32 as c_uint;
                j = 0i32;
                while j < 256i32 {
                    if 0 != *F_i_.offset(j as isize) {
                        if 0 != rle_j {
                            rle_j = rle_j.wrapping_sub(1)
                        } else {
                            let fresh49 = cp;
                            cp = cp.offset(1);
                            *fresh49 = j as c_uchar;
                            if 0 == rle_j && 0 != j &&
                                   0 != *F_i_.offset((j - 1i32) as isize) {
                                rle_j = (j + 1i32) as c_uint;
                                while rle_j < 256i32 as c_uint &&
                                          0 != *F_i_.offset(rle_j as isize) {
                                    rle_j = rle_j.wrapping_add(1)
                                }
                                rle_j =
                                    rle_j.wrapping_sub((j + 1i32) as
                                                           c_uint);
                                let fresh50 = cp;
                                cp = cp.offset(1);
                                *fresh50 = rle_j as c_uchar
                            }
                        }
                        if *F_i_.offset(j as isize) < 128i32 {
                            let fresh51 = cp;
                            cp = cp.offset(1);
                            *fresh51 =
                                *F_i_.offset(j as isize) as c_uchar
                        } else {
                            let fresh52 = cp;
                            cp = cp.offset(1);
                            *fresh52 =
                                (128i32 | *F_i_.offset(j as isize) >> 8i32) as
                                    c_uchar;
                            let fresh53 = cp;
                            cp = cp.offset(1);
                            *fresh53 =
                                (*F_i_.offset(j as isize) & 255i32) as
                                    c_uchar
                        }
                        Rans64EncSymbolInit(&mut syms[i as usize][j as usize]
                                            , x,
                                            *F_i_.offset(j as isize) as
                                                uint32_t, 12i32 as uint32_t);
                        x =
                            x.wrapping_add(*F_i_.offset(j as isize) as
                                               c_uint)
                    }
                    j += 1
                }
                let fresh54 = cp;
                cp = cp.offset(1);
                *fresh54 = 0i32 as c_uchar
            }
            i += 1
        }
        let fresh55 = cp;
        cp = cp.offset(1);
        *fresh55 = 0i32 as c_uchar;
        tab_size =
            out_buf.offset_to(cp).expect("bad offset_to") as c_long as
                c_uint;
        let mut rans: [RansState;4] = [0;4];
        Rans64EncInit(&mut rans[0]);
        Rans64EncInit(&mut rans[1]);
        Rans64EncInit(&mut rans[2]);
        Rans64EncInit(&mut rans[3]);
        let mut ptr: *mut uint8_t = out_end;
        let mut isz4: c_int = (in_size >> 2i32) as c_int;
        let mut i3: c_int = 4i32 * isz4 - 2i32;
        let mut l3: c_uchar = *in_0.offset((i3 + 1i32) as isize);
        let pptr = core::slice::from_raw_parts_mut(out_buf.offset(out_end_byte_offset as isize &3) as *mut uint32_t, out_end_byte_offset>>2);
        let mut ptr_offset = (out_end_byte_offset >> 2) - 1;
        l3 =
            *in_0.offset(in_size.wrapping_sub(1i32 as c_uint) as isize);
        i3 = in_size.wrapping_sub(2i32 as c_uint) as c_int;
        while i3 > 4i32 * isz4 - 2i32 {
            let mut c3: c_uchar = *in_0.offset(i3 as isize);
            //eprintln!("State {:x}",rans[3&THREE]);
            Rans64EncPutSymbol(&mut rans[3& THREE],
                               pptr, &mut ptr_offset,
                               &mut syms[c3 as usize][l3 as usize]
                                  , 12i32 as uint32_t);
            //eprintln!("PUT {:x} {:x} {:?}", c3, l3, syms[c3 as usize][l3 as usize]);
            //eprintln!("?tate {:x}",rans[3&THREE]);
            l3 = c3;
            i3 -= 1
        }
        assert_eq!(i3 & 3, 2);
        while i3 >= 0 {
            let mut c3: c_uchar = *in_0.offset(i3 as isize);
            //eprintln!("State {:x}[{}]",rans[(i3+1) as usize & THREE], i3 as usize & 3);
            Rans64EncPutSymbol(&mut rans[(i3+1) as usize & THREE],
                               pptr, &mut ptr_offset,
                               &mut syms[c3 as usize][l3 as usize]
                                  , 12i32 as uint32_t);
            //eprintln!("PUT {:x} {:x}  {:?}", c3, l3, syms[c3 as usize][l3 as usize]);
            //eprintln!("?tate {:x}[{}]",rans[(i3+1) as usize &THREE], i3 as usize & 3);
            l3 = c3;
            if i3 == 0 {
                break;
            }
            i3 = i3.wrapping_sub(1)
        }
        //eprintln!("State {:x}",rans[0&THREE]);
        Rans64EncPutSymbol(&mut rans[0&THREE],
                           pptr, &mut ptr_offset,
                           &mut syms[0usize][l3 as usize], 12i32 as uint32_t);
        //eprintln!("PUT {:x} {:x} {:?}", 0, l3, syms[0][l3 as usize]);
        //eprintln!("?tate {:x}",rans[0&THREE]);
        Rans64EncFlush(&mut rans[3],
                       pptr, &mut ptr_offset);
        Rans64EncFlush(&mut rans[2],
                       pptr, &mut ptr_offset);
        Rans64EncFlush(&mut rans[1],
                       pptr, &mut ptr_offset);
        Rans64EncFlush(&mut rans[0],
                       pptr, &mut ptr_offset);
        let out_byte_size = (pptr.len() - (ptr_offset.wrapping_add(1))) << 2;
        *out_size = out_byte_size as u32 + tab_size;
        cp = out_buf;
        let fresh56 = cp;
        cp = cp.offset(1);
        *fresh56 =
            (in_size >> 0i32 & 255i32 as c_uint) as c_uchar;
        let fresh57 = cp;
        cp = cp.offset(1);
        *fresh57 =
            (in_size >> 8i32 & 255i32 as c_uint) as c_uchar;
        let fresh58 = cp;
        cp = cp.offset(1);
        *fresh58 =
            (in_size >> 16i32 & 255i32 as c_uint) as c_uchar;
        let fresh59 = cp;
        cp = cp.offset(1);
        *fresh59 =
            (in_size >> 24i32 & 255i32 as c_uint) as c_uchar;
        memmove(out_buf.offset(tab_size as isize) as *mut c_void,
                pptr.as_ptr().offset(ptr_offset as isize + 1) as *const c_void,
                *out_size as c_ulong);
        return out_buf
    };
}



#[no_mangle]
pub unsafe extern "C" fn rans_uncompress_O1(mut in_0: *mut c_uchar,
                                            mut in_size: c_uint,
                                            mut out_size: *mut c_uint)
                                            -> *mut c_uchar {
    let mut cp: *mut c_uchar = in_0.offset(4isize);
    let mut ht_in_offset = 4;
    let mut i: c_int = 0;
    let mut j: c_int = -999i32;
    let mut x: c_int = 0;
    let mut out_sz: c_int = 0;
    let mut rle_i: c_int = 0;
    let mut rle_j: c_int = 0;
    let mut out_buf: *mut c_char = 0 as *mut c_char;
    let mut D: [ari_decoder; 256] = [ari_decoder{R: [0; 4096],}; 256];
    let mut syms: [[RansDecSymbol; 256]; 256] =
        [[RansDecSymbol{start: 0, freq: 0,}; 256]; 256];
    out_sz =
        (*in_0.offset(0isize) as c_int) << 0i32 |
            (*in_0.offset(1isize) as c_int) << 8i32 |
            (*in_0.offset(2isize) as c_int) << 16i32 |
            (*in_0.offset(3isize) as c_int) << 24i32;
    out_buf = malloc(out_sz as c_ulong) as *mut c_char;
    if !out_buf.is_null() {
        safe_rans_uncompress_O1(core::slice::from_raw_parts(in_0, in_size as usize),
                                core::slice::from_raw_parts_mut(out_buf, out_sz as usize));
    }
    *out_size = out_sz as c_uint;
    out_buf
}
fn safe_rans_uncompress_O1(in_0: &[c_uchar],
                           mut out_buf: &mut[c_uchar]) {
    let mut cp = &in_0[4..];
    let mut ht_in_offset = 4;
    let mut i: c_int = 0;
    let mut j: c_int = -999i32;
    let mut x: c_int = 0;
    let mut out_sz: c_int = 0;
    let mut rle_i: c_int = 0;
    let mut rle_j: c_int = 0;
    let mut D: [ari_decoder; 256] = [ari_decoder{R: [0; 4096],}; 256];
    let mut syms: [[RansDecSymbol; 256]; 256] =
        [[RansDecSymbol{start: 0, freq: 0,}; 256]; 256];
    out_sz =
        (in_0[0] as c_int) << 0i32 |
        (in_0[1] as c_int) << 8i32 |
        (in_0[2] as c_int) << 16i32 |
        (in_0[3] as c_int) << 24i32;
    assert_eq!(out_sz as usize, out_buf.len());
    rle_i = 0i32;
    let fresh60 = cp;
    cp = &cp[1..];
    ht_in_offset += 1;
    i = fresh60[0] as c_int;
    loop {
        x = 0i32;
        rle_j = x;
        let fresh61 = cp;
        cp = &cp[1..];
        ht_in_offset += 1;
        j = fresh61[0] as c_int;
        loop {
            let mut F: c_int = 0;
            let mut C: c_int = 0;
            let fresh62 = cp;
            cp = &cp[1..];
            ht_in_offset += 1;
            F = fresh62[0] as c_int;
            if F >= 128i32 {
                F &= !128i32;
                let fresh63 = cp;
                cp = &cp[1..];
                ht_in_offset += 1;
                F = (F & 127i32) << 8i32 | fresh63[0] as c_int
            }
            C = x;
            if 0 == F { F = 1i32 << 12i32 }
            Rans64DecSymbolInit(&mut syms[i as usize][j as usize]
                                , C as u16,
                                F as u16);
            for item in D[i as usize].R.split_at_mut(x as usize).1.split_at_mut(F as usize).0.iter_mut() {
                *item = j as u8;
            }
            //memset(&mut  as *mut c_uchar
            //       as *mut c_void, j, F as c_ulong);
            x += F;
            if 0 == rle_j && j + 1i32 == cp[0] as c_int {
                let fresh64 = cp;
                cp = &cp[1..];
                ht_in_offset += 1;
                j = fresh64[0] as c_int;
                let fresh65 = cp;
                cp = &cp[1..];
                ht_in_offset += 1;
                rle_j = fresh65[0] as c_int
            } else if 0 != rle_j {
                rle_j -= 1;
                j += 1
            } else {
                let fresh66 = cp;
                cp = &cp[1..];
                ht_in_offset += 1;
                j = fresh66[0] as c_int
            }
            if !(0 != j) { break ; }
        }
        if 0 == rle_i && i + 1i32 == cp[0] as c_int {
            let fresh67 = cp;
            cp = &cp[1..];
            ht_in_offset += 1;
            i = fresh67[0] as c_int;
            let fresh68 = cp;
            cp = &cp[1..];
            ht_in_offset += 1;
            rle_i = fresh68[0] as c_int
        } else if 0 != rle_i {
            rle_i -= 1;
            i += 1
        } else {
            let fresh69 = cp;
            cp = &cp[1..];
            ht_in_offset += 1;
            i = fresh69[0] as i32;
        }
        if !(0 != i) { break ; }
    }
    let mut rans0: RansState = 0;
    let mut rans1: RansState = 0;
    let mut rans2: RansState = 0;
    let mut rans3: RansState = 0;

    let pptr = cp;
    let mut ptr_offset = 0u32;
    u8Rans64DecInit(&mut rans0,
                  pptr, &mut ptr_offset);
    
    u8Rans64DecInit(&mut rans1,
                  pptr, &mut ptr_offset);
    u8Rans64DecInit(&mut rans2,
                  pptr, &mut ptr_offset);
    u8Rans64DecInit(&mut rans3,
                  pptr, &mut ptr_offset);
    let mut isz4: c_int = out_sz >> 2i32;
    let mut l0: c_int = 0i32;
    let mut l1: c_int = 0i32;
    let mut l2: c_int = 0i32;
    let mut l3: c_int = 0i32;
    let mut i4: [c_int; 4] =
        [0i32 * isz4, 1i32 * isz4, 2i32 * isz4, 3i32 * isz4];
    let mut R: [RansState; 4] = [0; 4];
    R[0usize] = rans0;
    R[1usize] = rans1;
    R[2usize] = rans2;
    R[3usize] = rans3;
    {
    let (out_buf0, outrest0) = out_buf.split_at_mut(isz4 as usize);
    let (out_buf1, outrest1) = outrest0.split_at_mut(isz4 as usize);
    let (out_buf2, out_buf3) = outrest1.split_at_mut(isz4 as usize);
    for (out0, (out1, (out2, out3))) in out_buf0.iter_mut().zip(out_buf1.iter_mut().zip(out_buf2.iter_mut().zip(out_buf3.iter_mut()))) {
        let mut m: [uint32_t; 4] =
            [(R[0usize] &
              (1u32 << 12i32).wrapping_sub(1i32 as c_uint) as
              c_ulong) as uint32_t,
             (R[1usize] &
              (1u32 << 12i32).wrapping_sub(1i32 as c_uint) as
              c_ulong) as uint32_t,
             (R[2usize] &
              (1u32 << 12i32).wrapping_sub(1i32 as c_uint) as
              c_ulong) as uint32_t,
             (R[3usize] &
              (1u32 << 12i32).wrapping_sub(1i32 as c_uint) as
              c_ulong) as uint32_t];
        let mut c: [uint8_t; 4] =
            [D[l0 as usize].R[m[0usize] as usize],
             D[l1 as usize].R[m[1usize] as usize],
             D[l2 as usize].R[m[2usize] as usize],
             D[l3 as usize].R[m[3usize] as usize]];
        *out0 = c[0usize] as c_char;
        *out1 = c[1usize] as c_char;
        *out2 = c[2usize] as c_char;
        *out3 = c[3usize] as c_char;
        R[0usize] =
            (syms[l0 as usize][c[0usize] as usize].freq as
             c_ulong).wrapping_mul(R[0usize] >> 12i32);
        R[1usize] =
            (syms[l1 as usize][c[1usize] as usize].freq as
             c_ulong).wrapping_mul(R[1usize] >> 12i32);
        R[2usize] =
            (syms[l2 as usize][c[2usize] as usize].freq as
             c_ulong).wrapping_mul(R[2usize] >> 12i32);
        R[3usize] =
            (syms[l3 as usize][c[3usize] as usize].freq as
             c_ulong).wrapping_mul(R[3usize] >> 12i32);
        R[0usize] =
            (R[0usize] as
             c_ulong).wrapping_add(m[0usize].wrapping_sub(u32::from(syms[l0
                                                               as
                                                               usize][c[0usize]
                                                                      as
                                                                      usize].start))
                                   as c_ulong) as
            RansState as RansState;
        R[1usize] =
            (R[1usize] as
             c_ulong).wrapping_add(m[1usize].wrapping_sub(u32::from(syms[l1
                                                               as
                                                               usize][c[1usize]
                                                                      as
                                                                      usize].start))
                                   as c_ulong) as
            RansState as RansState;
        R[2usize] =
            (R[2usize] as
             c_ulong).wrapping_add(m[2usize].wrapping_sub(u32::from(syms[l2
                                                                                 as
                                                               usize][c[2usize]
                                                                      as
                                                                                            usize].start))
                                   as c_ulong) as
            RansState as RansState;
        R[3usize] =
            (R[3usize] as
             c_ulong).wrapping_add(m[3usize].wrapping_sub(u32::from(syms[l3
                                                               as
                                                               usize][c[3usize]
                                                                      as
                                                                      usize].start))
                                   as c_ulong) as
            RansState as RansState;
        u8Rans64DecRenorm(&mut R[0usize],
                          pptr, &mut ptr_offset);
        u8Rans64DecRenorm(&mut R[1usize],
                          pptr, &mut ptr_offset);
        u8Rans64DecRenorm(&mut R[2usize],
                          pptr, &mut ptr_offset);
        u8Rans64DecRenorm(&mut R[3usize],
                          pptr, &mut ptr_offset);
        l0 = c[0usize] as c_int;
        l1 = c[1usize] as c_int;
        l2 = c[2usize] as c_int;
        l3 = c[3usize] as c_int;
        i4[0usize] += 1;
        i4[1usize] += 1;
        i4[2usize] += 1;
        i4[3usize] += 1
    }
    }
    rans0 = R[0usize];
    rans1 = R[1usize];
    rans2 = R[2usize];
    rans3 = R[3usize];
    while i4[3usize] < out_sz {
        let mut c3: c_uchar =
            D[l3 as
              usize].R[Rans64DecGet(&mut rans3,
                                    12i32 as uint32_t) as usize];
        out_buf[i4[3usize] as usize] = c3 as c_char;
        u8Rans64DecAdvanceSymbol(&mut rans3,
                                 pptr, &mut ptr_offset,
                                 &mut syms[l3 as usize][c3 as usize], 12i32 as uint32_t);
        l3 = c3 as c_int;
        i4[3usize] += 1
    }
    //*out_size = out_sz as c_uint;
}












#[no_mangle]
pub unsafe extern "C" fn rans_uncompress_O1d(mut in_0: *mut c_uchar,
                                            mut in_size: c_uint,
                                            mut out_size: *mut c_uint)
                                            -> *mut c_uchar {
    let mut cp: *mut c_uchar = in_0.offset(4isize);
    let mut ht_in_offset = 4;
    let mut i: c_int = 0;
    let mut j: c_int = -999i32;
    let mut x: c_int = 0;
    let mut out_sz: c_int = 0;
    let mut rle_i: c_int = 0;
    let mut rle_j: c_int = 0;
    let mut out_buf: *mut c_char = 0 as *mut c_char;
    let mut D: [ari_decoder; 256] = [ari_decoder{R: [0; 4096],}; 256];
    let mut syms: [[RansDecSymbol; 256]; 256] =
        [[RansDecSymbol{start: 0, freq: 0,}; 256]; 256];
    out_sz =
        (*in_0.offset(0isize) as c_int) << 0i32 |
            (*in_0.offset(1isize) as c_int) << 8i32 |
            (*in_0.offset(2isize) as c_int) << 16i32 |
            (*in_0.offset(3isize) as c_int) << 24i32;
    out_buf = malloc(out_sz as c_ulong) as *mut c_char;
    if !out_buf.is_null() {
        safe_rans_uncompress_O1d(core::slice::from_raw_parts(in_0, in_size as usize),
                                core::slice::from_raw_parts_mut(out_buf, out_sz as usize));
    }
    *out_size = out_sz as c_uint;
    out_buf
}
fn safe_rans_uncompress_O1d(in_0: &[c_uchar],
                            mut out_buf: &mut[c_uchar]) {
    let mut cp = &in_0[4..];
    let mut ht_in_offset = 4;
    let mut i: c_int = 0;
    let mut j: c_int = -999i32;
    let mut x: c_int = 0;
    let mut out_sz: c_int = 0;
    let mut rle_i: c_int = 0;
    let mut rle_j: c_int = 0;
    let mut D: [ari_decoder; 256] = [ari_decoder{R: [0; 4096],}; 256];
    let mut syms: [[RansDecSymbol; 256]; 256] =
        [[RansDecSymbol{start: 0, freq: 0,}; 256]; 256];
    out_sz =
        (in_0[0] as c_int) << 0i32 |
        (in_0[1] as c_int) << 8i32 |
        (in_0[2] as c_int) << 16i32 |
        (in_0[3] as c_int) << 24i32;
    assert_eq!(out_sz as usize, out_buf.len());
    rle_i = 0i32;
    let fresh60 = cp;
    cp = &cp[1..];
    ht_in_offset += 1;
    i = fresh60[0] as c_int;
    loop {
        x = 0i32;
        rle_j = x;
        let fresh61 = cp;
        cp = &cp[1..];
        ht_in_offset += 1;
        j = fresh61[0] as c_int;
        loop {
            let mut F: c_int = 0;
            let mut C: c_int = 0;
            let fresh62 = cp;
            cp = &cp[1..];
            ht_in_offset += 1;
            F = fresh62[0] as c_int;
            if F >= 128i32 {
                F &= !128i32;
                let fresh63 = cp;
                cp = &cp[1..];
                ht_in_offset += 1;
                F = (F & 127i32) << 8i32 | fresh63[0] as c_int
            }
            C = x;
            if 0 == F { F = 1i32 << 12i32 }
            Rans64DecSymbolInit(&mut syms[i as usize][j as usize]
                                , C as u16,
                                F as u16);
            for item in D[i as usize].R.split_at_mut(x as usize).1.split_at_mut(F as usize).0.iter_mut() {
                *item = j as u8;
            }
            //memset(&mut  as *mut c_uchar
            //       as *mut c_void, j, F as c_ulong);
            x += F;
            if 0 == rle_j && j + 1i32 == cp[0] as c_int {
                let fresh64 = cp;
                cp = &cp[1..];
                ht_in_offset += 1;
                j = fresh64[0] as c_int;
                let fresh65 = cp;
                cp = &cp[1..];
                ht_in_offset += 1;
                rle_j = fresh65[0] as c_int
            } else if 0 != rle_j {
                rle_j -= 1;
                j += 1
            } else {
                let fresh66 = cp;
                cp = &cp[1..];
                ht_in_offset += 1;
                j = fresh66[0] as c_int
            }
            if !(0 != j) { break ; }
        }
        if 0 == rle_i && i + 1i32 == cp[0] as c_int {
            let fresh67 = cp;
            cp = &cp[1..];
            ht_in_offset += 1;
            i = fresh67[0] as c_int;
            let fresh68 = cp;
            cp = &cp[1..];
            ht_in_offset += 1;
            rle_i = fresh68[0] as c_int
        } else if 0 != rle_i {
            rle_i -= 1;
            i += 1
        } else {
            let fresh69 = cp;
            cp = &cp[1..];
            ht_in_offset += 1;
            i = fresh69[0] as i32;
        }
        if !(0 != i) { break ; }
    }
    let mut R: [RansState; 4] = [0; 4];

    let pptr = cp;
    let mut ptr_offset = 0u32;
    u8Rans64DecInit(&mut R[0],
                  pptr, &mut ptr_offset);
    
    u8Rans64DecInit(&mut R[1],
                  pptr, &mut ptr_offset);
    u8Rans64DecInit(&mut R[2],
                  pptr, &mut ptr_offset);
    u8Rans64DecInit(&mut R[3],
                    pptr, &mut ptr_offset);
    ////eprintln!("{:x}\n{:x}\n{:x}\n{:x}", R[0], R[1], R[2], R[3]);
    let mut isz4: c_int = out_sz >> 2i32;
    let mut l0: c_int = 0i32;
    let mut i4: [c_int; 1] =
        [0i32];
    const prob_mask:u16 = (1u16 << 12) - 1;
    let mut m = R[0] as u16 & prob_mask;
    let mut next_fetch = &D[l0 as usize].R[m as usize];
    let mut next_c = *next_fetch;
    for (index, out0) in out_buf.split_at_mut(isz4 as usize * 4).0.iter_mut().enumerate() {
        //eprintln!("State {:x}[{}]", R[index & THREE], index & THREE);
        let mut c = next_c;
        *out0 = c;
        let start_freq = syms[l0 as usize][c as usize];
        let mnew = R[(index + 1)& THREE] as u16 & prob_mask;
        R[index & THREE] =
            (start_freq.freq as
             c_ulong).wrapping_mul(R[index & THREE] >> 12);
        R[index & THREE] =
            (R[index & THREE] as
             c_ulong).wrapping_add(m.wrapping_sub(start_freq.start) as u64);
        let lt_1_sl_31 = (R[index & THREE] & 0xffff_ffff_8000_0000);
        next_fetch = &D[c as usize].R[mnew as usize];
        next_c = *next_fetch;
        m = mnew;
        if lt_1_sl_31 == 0 {
            u8Rans64DecForceRenorm(&mut R[index & THREE],
                              pptr, &mut ptr_offset);
        }
        //eprintln!("GET {:x} {:x} {:?}", l0, c, syms[l0 as usize][c as usize]);
        //eprintln!("Atate {:x}[{}]", R[index & THREE], index & THREE);
        l0 = c as c_int;
        i4[0usize] += 1;
    }
    while i4[0usize] < out_sz {
        //eprintln!("State {:x}", R[3 & THREE]);
        let mut c3: c_uchar =
            D[l0 as
              usize].R[Rans64DecGet(&mut R[3 & THREE],
                                    12 as uint32_t) as usize];
        out_buf[i4[0usize] as usize] = c3 as c_char;
        u8Rans64DecAdvanceSymbol(&mut R[THREE],
                                 pptr, &mut ptr_offset,
                                 &mut syms[l0 as usize][c3 as usize], 12i32 as uint32_t);
        //eprintln!("GET {:x} {:x}", l0, c3);
        //eprintln!("Atate {:x}", R[THREE]);
        l0 = c3 as c_int;
        i4[0usize] += 1
    }
    //*out_size = out_sz as c_uint;
}




fn safe_rans_uncompress_O1g(in_0: &[c_uchar],
                            mut out_buf: &mut[c_uchar]) {
    let mut cp = &in_0[4..];
    let mut ht_in_offset = 4;
    let mut i: c_int = 0;
    let mut j: c_int = -999i32;
    let mut x: c_int = 0;
    let mut out_sz: c_int = 0;
    let mut rle_i: c_int = 0;
    let mut rle_j: c_int = 0;
    let mut D: [ari_decoder; 256] = [ari_decoder{R: [0; 4096],}; 256];
    let mut syms: [[RansDecSymbol; 256]; 256] =
        [[RansDecSymbol{start: 0, freq: 0,}; 256]; 256];
    out_sz =
        (in_0[0] as c_int) << 0i32 |
        (in_0[1] as c_int) << 8i32 |
        (in_0[2] as c_int) << 16i32 |
        (in_0[3] as c_int) << 24i32;
    assert_eq!(out_sz as usize, out_buf.len());
    rle_i = 0i32;
    let fresh60 = cp;
    cp = &cp[1..];
    ht_in_offset += 1;
    i = fresh60[0] as c_int;
    loop {
        x = 0i32;
        rle_j = x;
        let fresh61 = cp;
        cp = &cp[1..];
        ht_in_offset += 1;
        j = fresh61[0] as c_int;
        loop {
            let mut F: c_int = 0;
            let mut C: c_int = 0;
            let fresh62 = cp;
            cp = &cp[1..];
            ht_in_offset += 1;
            F = fresh62[0] as c_int;
            if F >= 128i32 {
                F &= !128i32;
                let fresh63 = cp;
                cp = &cp[1..];
                ht_in_offset += 1;
                F = (F & 127i32) << 8i32 | fresh63[0] as c_int
            }
            C = x;
            if 0 == F { F = 1i32 << 12i32 }
            Rans64DecSymbolInit(&mut syms[i as usize][j as usize]
                                , C as u16,
                                F as u16);
            for item in D[i as usize].R.split_at_mut(x as usize).1.split_at_mut(F as usize).0.iter_mut() {
                *item = j as u8;
            }
            //memset(&mut  as *mut c_uchar
            //       as *mut c_void, j, F as c_ulong);
            x += F;
            if 0 == rle_j && j + 1i32 == cp[0] as c_int {
                let fresh64 = cp;
                cp = &cp[1..];
                ht_in_offset += 1;
                j = fresh64[0] as c_int;
                let fresh65 = cp;
                cp = &cp[1..];
                ht_in_offset += 1;
                rle_j = fresh65[0] as c_int
            } else if 0 != rle_j {
                rle_j -= 1;
                j += 1
            } else {
                let fresh66 = cp;
                cp = &cp[1..];
                ht_in_offset += 1;
                j = fresh66[0] as c_int
            }
            if !(0 != j) { break ; }
        }
        if 0 == rle_i && i + 1i32 == cp[0] as c_int {
            let fresh67 = cp;
            cp = &cp[1..];
            ht_in_offset += 1;
            i = fresh67[0] as c_int;
            let fresh68 = cp;
            cp = &cp[1..];
            ht_in_offset += 1;
            rle_i = fresh68[0] as c_int
        } else if 0 != rle_i {
            rle_i -= 1;
            i += 1
        } else {
            let fresh69 = cp;
            cp = &cp[1..];
            ht_in_offset += 1;
            i = fresh69[0] as i32;
        }
        if !(0 != i) { break ; }
    }
    let mut next_state = [[RansDecSymbolStartFreq{start:0,freq:0,sym:0};4096];256];
    for prior in 0..256 {
        for state in 0..4096 {
            let sym = D[prior].R[state];
            next_state[prior][state] = RansDecSymbolStartFreq{
                sym: sym,
                start:syms[prior][sym as usize].start as u16,
                freq:syms[prior][sym as usize].freq as u16,
            };
        }
    }
    let mut R: [RansState; 4] = [0; 4];

    let pptr = cp;
    let mut ptr_offset = 0u32;
    u8Rans64DecInit(&mut R[0],
                  pptr, &mut ptr_offset);
    
    u8Rans64DecInit(&mut R[1],
                  pptr, &mut ptr_offset);
    u8Rans64DecInit(&mut R[2],
                  pptr, &mut ptr_offset);
    u8Rans64DecInit(&mut R[3],
                    pptr, &mut ptr_offset);
    ////eprintln!("{:x}\n{:x}\n{:x}\n{:x}", R[0], R[1], R[2], R[3]);
    let mut isz4: c_int = out_sz >> 2i32;
    let mut l0: c_int = 0i32;
    let mut i4: [c_int; 1] =
        [0i32];
    const prob_mask:u16 = (1u16 << 12) - 1;
    for (index, out0) in out_buf.split_at_mut(isz4 as usize * 4).0.iter_mut().enumerate() {
        //eprintln!("State {:x}[{}]", R[index & THREE], index & THREE);
        let mut m = R[index & THREE] as u16& prob_mask;
        let mut c_start_freq = next_state[l0 as usize][m as usize];
        *out0 = c_start_freq.sym;
        R[index & THREE] =
            (c_start_freq.freq as
             c_ulong).wrapping_mul(R[index & THREE] >> 12);
        R[index & THREE] =
            (R[index & THREE] as
             c_ulong).wrapping_add(m.wrapping_sub(c_start_freq.start)
                                   as c_ulong) as
            RansState as RansState;
        u8Rans64DecRenorm(&mut R[index & THREE],
                          pptr, &mut ptr_offset);
        //eprintln!("GET {:x} {:x} {:?}", l0, c, syms[l0 as usize][c as usize]);
        //eprintln!("Atate {:x}[{}]", R[index & THREE], index & THREE);
        l0 = c_start_freq.sym as c_int;
        i4[0usize] += 1;
    }
    while i4[0usize] < out_sz {
        //eprintln!("State {:x}", R[3 & THREE]);
        let mut c3: c_uchar =
            D[l0 as
              usize].R[Rans64DecGet(&mut R[3 & THREE],
                                    12 as uint32_t) as usize];
        out_buf[i4[0usize] as usize] = c3 as c_char;
        u8Rans64DecAdvanceSymbol(&mut R[THREE],
                                 pptr, &mut ptr_offset,
                                 &mut syms[l0 as usize][c3 as usize], 12i32 as uint32_t);
        //eprintln!("GET {:x} {:x}", l0, c3);
        //eprintln!("Atate {:x}", R[THREE]);
        l0 = c3 as c_int;
        i4[0usize] += 1
    }
    //*out_size = out_sz as c_uint;
}

fn safe_rans_uncompress_O1e(in_0: &[c_uchar],
                            mut out_buf: &mut[c_uchar]) {
    assert_eq!(THREE, 3);
    let mut cp = &in_0[4..];
    let mut ht_in_offset = 4;
    let mut i: c_int = 0;
    let mut j: c_int = -999i32;
    let mut x: c_int = 0;
    let mut out_sz: c_int = 0;
    let mut rle_i: c_int = 0;
    let mut rle_j: c_int = 0;
    let mut D: [ari_decoder; 256] = [ari_decoder{R: [0; 4096],}; 256];
    let mut syms: [[RansDecSymbol; 256]; 256] =
        [[RansDecSymbol{start: 0, freq: 0,}; 256]; 256];
    out_sz =
        (in_0[0] as c_int) << 0i32 |
        (in_0[1] as c_int) << 8i32 |
        (in_0[2] as c_int) << 16i32 |
        (in_0[3] as c_int) << 24i32;
    assert_eq!(out_sz as usize, out_buf.len());
    rle_i = 0i32;
    let fresh60 = cp;
    cp = &cp[1..];
    ht_in_offset += 1;
    i = fresh60[0] as c_int;
    loop {
        x = 0i32;
        rle_j = x;
        let fresh61 = cp;
        cp = &cp[1..];
        ht_in_offset += 1;
        j = fresh61[0] as c_int;
        loop {
            let mut F: c_int = 0;
            let mut C: c_int = 0;
            let fresh62 = cp;
            cp = &cp[1..];
            ht_in_offset += 1;
            F = fresh62[0] as c_int;
            if F >= 128i32 {
                F &= !128i32;
                let fresh63 = cp;
                cp = &cp[1..];
                ht_in_offset += 1;
                F = (F & 127i32) << 8i32 | fresh63[0] as c_int
            }
            C = x;
            if 0 == F { F = 1i32 << 12i32 }
            Rans64DecSymbolInit(&mut syms[i as usize][j as usize]
                                , C as u16,
                                F as u16);
            for item in D[i as usize].R.split_at_mut(x as usize).1.split_at_mut(F as usize).0.iter_mut() {
                *item = j as u8;
            }
            //memset(&mut  as *mut c_uchar
            //       as *mut c_void, j, F as c_ulong);
            x += F;
            if 0 == rle_j && j + 1i32 == cp[0] as c_int {
                let fresh64 = cp;
                cp = &cp[1..];
                ht_in_offset += 1;
                j = fresh64[0] as c_int;
                let fresh65 = cp;
                cp = &cp[1..];
                ht_in_offset += 1;
                rle_j = fresh65[0] as c_int
            } else if 0 != rle_j {
                rle_j -= 1;
                j += 1
            } else {
                let fresh66 = cp;
                cp = &cp[1..];
                ht_in_offset += 1;
                j = fresh66[0] as c_int
            }
            if !(0 != j) { break ; }
        }
        if 0 == rle_i && i + 1i32 == cp[0] as c_int {
            let fresh67 = cp;
            cp = &cp[1..];
            ht_in_offset += 1;
            i = fresh67[0] as c_int;
            let fresh68 = cp;
            cp = &cp[1..];
            ht_in_offset += 1;
            rle_i = fresh68[0] as c_int
        } else if 0 != rle_i {
            rle_i -= 1;
            i += 1
        } else {
            let fresh69 = cp;
            cp = &cp[1..];
            ht_in_offset += 1;
            i = fresh69[0] as i32;
        }
        if !(0 != i) { break ; }
    }
    let mut R: [RansState; 4] = [0; 4];

    let pptr = cp;
    let mut ptr_offset = 0u32;
    u8Rans64DecInit(&mut R[0],
                  pptr, &mut ptr_offset);
    
    u8Rans64DecInit(&mut R[1],
                  pptr, &mut ptr_offset);
    u8Rans64DecInit(&mut R[2],
                  pptr, &mut ptr_offset);
    u8Rans64DecInit(&mut R[3],
                    pptr, &mut ptr_offset);
    ////eprintln!("{:x}\n{:x}\n{:x}\n{:x}", R[0], R[1], R[2], R[3]);
    let mut isz4: c_int = out_sz >> 2i32;
    let mut l0: c_int = 0i32;
    const prob_mask:u16 = (1u16 << 12) - 1;
    let mut last_d = &D[l0 as usize];
    assert!(isz4 as usize * 4 <= out_buf.len());
    for index in 0..isz4 {
        assert!(index as usize * 4 + 3 <= out_buf.len());
        let mut m = [R[0] as u16 & prob_mask,
                     R[1] as u16 & prob_mask,
                     R[2] as u16 & prob_mask,
                     R[3] as u16 & prob_mask];
        let mut c0 = last_d.R[m[0] as usize];
        
        let mut c1 = D[c0 as usize].R[m[1] as usize];
        let mut c2 = D[c1 as usize].R[m[2] as usize];
        let mut c3 = D[c2 as usize].R[m[3] as usize];
        //eprintln!("State {:x}[{}]", R[0], 0);
        let sym = [syms[l0 as usize][c0 as usize],
                   syms[c0 as usize][c1 as usize],
                   syms[c1 as usize][c2 as usize],
                   syms[c2 as usize][c3 as usize]];
        last_d = &D[c3 as usize];
        R[0] = (sym[0].freq as
             c_ulong).wrapping_mul(R[0] >> 12);
        R[1] =
            (sym[1].freq as
             c_ulong).wrapping_mul(R[1] >> 12);
        R[2] =
            (sym[2].freq as
             c_ulong).wrapping_mul(R[2] >> 12);
        R[3] =
            (sym[3].freq as
             c_ulong).wrapping_mul(R[3] >> 12);
        R[0] = R[0].wrapping_add(m[0].wrapping_sub(sym[0].start as u16)
                                 as c_ulong);
        R[1] = R[1].wrapping_add(m[1].wrapping_sub(sym[1].start as u16)
                                 as c_ulong);
        R[2] = R[2].wrapping_add(m[2].wrapping_sub(sym[2].start as u16)
                                 as c_ulong);
        R[3] = R[3].wrapping_add(m[3].wrapping_sub(sym[3].start as u16)
                                 as c_ulong);
        u8Rans64DecRenorm(&mut R[0],
                          pptr, &mut ptr_offset);
        u8Rans64DecRenorm(&mut R[1],
                          pptr, &mut ptr_offset);
        u8Rans64DecRenorm(&mut R[2],
                          pptr, &mut ptr_offset);
        u8Rans64DecRenorm(&mut R[3],
                          pptr, &mut ptr_offset);
        //eprintln!("GET {:x} {:x} {:?}", l0, c0, sym[0]);
        //eprintln!("Atate {:x}[{}]", R[0], 0);
        //eprintln!("State {:x}[{}]", R[1], 1);
        //eprintln!("GET {:x} {:x} {:?}", c0, c1, sym[1]);
        //eprintln!("Atate {:x}[{}]", R[1], 1);
        //eprintln!("State {:x}[{}]", R[2], 2);
        //eprintln!("GET {:x} {:x} {:?}", c1, c2, sym[2]);
        //eprintln!("Atate {:x}[{}]", R[2], 2);
        //eprintln!("State {:x}[{}]", R[3], 3);
        //eprintln!("GET {:x} {:x} {:?}", c2, c3, sym[3]);
        //eprintln!("Atate {:x}[{}]", R[3], 3);
        out_buf[index as usize * 4] = c0;
        out_buf[index as usize * 4 + 1] = c1;
        out_buf[index as usize * 4 + 2] = c2;
        out_buf[index as usize * 4 + 3] = c3;
        l0 = c3 as c_int;
    }
    let mut i4 = isz4 * 4;
    while i4 < out_sz {
        //eprintln!("State {:x}", R[3 & THREE]);
        let mut c3: c_uchar =
            D[l0 as
              usize].R[Rans64DecGet(&mut R[3 & THREE],
                                    12 as uint32_t) as usize];
        out_buf[i4 as usize] = c3 as c_char;
        u8Rans64DecAdvanceSymbol(&mut R[THREE],
                                 pptr, &mut ptr_offset,
                                 &mut syms[l0 as usize][c3 as usize], 12i32 as uint32_t);
        //eprintln!("GET {:x} {:x}", l0, c3);
        //eprintln!("Atate {:x}", R[THREE]);
        l0 = c3 as c_int;
        i4 += 1
    }
    //*out_size = out_sz as c_uint;
}


fn safe_rans_uncompress_O1f(in_0: &[c_uchar],
                            mut out_buf: &mut[c_uchar]) {
    assert_eq!(THREE, 1);
    let mut cp = &in_0[4..];
    let mut ht_in_offset = 4;
    let mut i: c_int = 0;
    let mut j: c_int = -999i32;
    let mut x: c_int = 0;
    let mut out_sz: c_int = 0;
    let mut rle_i: c_int = 0;
    let mut rle_j: c_int = 0;
    let mut D: [ari_decoder; 256] = [ari_decoder{R: [0; 4096],}; 256];
    let mut syms: [[RansDecSymbol; 256]; 256] =
        [[RansDecSymbol{start: 0, freq: 0,}; 256]; 256];
    out_sz =
        (in_0[0] as c_int) << 0i32 |
        (in_0[1] as c_int) << 8i32 |
        (in_0[2] as c_int) << 16i32 |
        (in_0[3] as c_int) << 24i32;
    assert_eq!(out_sz as usize, out_buf.len());
    rle_i = 0i32;
    let fresh60 = cp;
    cp = &cp[1..];
    ht_in_offset += 1;
    i = fresh60[0] as c_int;
    loop {
        x = 0i32;
        rle_j = x;
        let fresh61 = cp;
        cp = &cp[1..];
        ht_in_offset += 1;
        j = fresh61[0] as c_int;
        loop {
            let mut F: c_int = 0;
            let mut C: c_int = 0;
            let fresh62 = cp;
            cp = &cp[1..];
            ht_in_offset += 1;
            F = fresh62[0] as c_int;
            if F >= 128i32 {
                F &= !128i32;
                let fresh63 = cp;
                cp = &cp[1..];
                ht_in_offset += 1;
                F = (F & 127i32) << 8i32 | fresh63[0] as c_int
            }
            C = x;
            if 0 == F { F = 1i32 << 12i32 }
            Rans64DecSymbolInit(&mut syms[i as usize][j as usize]
                                , C as u16,
                                F as u16);
            for item in D[i as usize].R.split_at_mut(x as usize).1.split_at_mut(F as usize).0.iter_mut() {
                *item = j as u8;
            }
            //memset(&mut  as *mut c_uchar
            //       as *mut c_void, j, F as c_ulong);
            x += F;
            if 0 == rle_j && j + 1i32 == cp[0] as c_int {
                let fresh64 = cp;
                cp = &cp[1..];
                ht_in_offset += 1;
                j = fresh64[0] as c_int;
                let fresh65 = cp;
                cp = &cp[1..];
                ht_in_offset += 1;
                rle_j = fresh65[0] as c_int
            } else if 0 != rle_j {
                rle_j -= 1;
                j += 1
            } else {
                let fresh66 = cp;
                cp = &cp[1..];
                ht_in_offset += 1;
                j = fresh66[0] as c_int
            }
            if !(0 != j) { break ; }
        }
        if 0 == rle_i && i + 1i32 == cp[0] as c_int {
            let fresh67 = cp;
            cp = &cp[1..];
            ht_in_offset += 1;
            i = fresh67[0] as c_int;
            let fresh68 = cp;
            cp = &cp[1..];
            ht_in_offset += 1;
            rle_i = fresh68[0] as c_int
        } else if 0 != rle_i {
            rle_i -= 1;
            i += 1
        } else {
            let fresh69 = cp;
            cp = &cp[1..];
            ht_in_offset += 1;
            i = fresh69[0] as i32;
        }
        if !(0 != i) { break ; }
    }
    let mut R: [RansState; 4] = [0; 4];

    let pptr = cp;
    let mut ptr_offset = 0u32;
    u8Rans64DecInit(&mut R[0],
                  pptr, &mut ptr_offset);
    
    u8Rans64DecInit(&mut R[1],
                  pptr, &mut ptr_offset);
    u8Rans64DecInit(&mut R[2],
                  pptr, &mut ptr_offset);
    u8Rans64DecInit(&mut R[3],
                    pptr, &mut ptr_offset);
    ////eprintln!("{:x}\n{:x}\n{:x}\n{:x}", R[0], R[1], R[2], R[3]);
    let mut isz4: c_int = out_sz >> 2i32;
    let mut l0: c_int = 0i32;
    const prob_mask:u16 = (1u16 << 12) - 1;
    let mut last_d = &D[l0 as usize];
    assert!(isz4 as usize * 4 <= out_buf.len());
    for index in 0..isz4*2 {
        assert!(index as usize * 2 + 1 <= out_buf.len());
        let mut m = [R[0] as u16 & prob_mask,
                     R[1] as u16 & prob_mask];
        let mut c0 = last_d.R[m[0] as usize];
        
        let mut c1 = D[c0 as usize].R[m[1] as usize];
        //eprintln!("State {:x}[{}]", R[0], 0);
        let sym = [syms[l0 as usize][c0 as usize],
                   syms[c0 as usize][c1 as usize]];
        last_d = &D[c1 as usize];
        R[0] = (sym[0].freq as
             c_ulong).wrapping_mul(R[0] >> 12);
        R[1] =
            (sym[1].freq as
             c_ulong).wrapping_mul(R[1] >> 12);
        R[0] = R[0].wrapping_add(m[0].wrapping_sub(sym[0].start as u16)
                                 as c_ulong);
        R[1] = R[1].wrapping_add(m[1].wrapping_sub(sym[1].start as u16)
                                 as c_ulong);
        u8Rans64DecRenorm(&mut R[0],
                          pptr, &mut ptr_offset);
        u8Rans64DecRenorm(&mut R[1],
                          pptr, &mut ptr_offset);
        //eprintln!("GET {:x} {:x} {:?}", l0, c0, sym[0]);
        //eprintln!("Atate {:x}[{}]", R[0], 0);
        //eprintln!("State {:x}[{}]", R[1], 1);
        //eprintln!("GET {:x} {:x} {:?}", c0, c1, sym[1]);
        //eprintln!("Atate {:x}[{}]", R[1], 1);
        //eprintln!("State {:x}[{}]", R[2], 2);
        //eprintln!("GET {:x} {:x} {:?}", c1, c2, sym[2]);
        //eprintln!("Atate {:x}[{}]", R[2], 2);
        //eprintln!("State {:x}[{}]", R[3], 3);
        //eprintln!("GET {:x} {:x} {:?}", c2, c3, sym[3]);
        //eprintln!("Atate {:x}[{}]", R[3], 3);
        out_buf[index as usize * 2] = c0;
        out_buf[index as usize * 2 + 1] = c1;
        l0 = c1 as c_int;
    }
    let mut i4 = isz4 * 4;
    while i4 < out_sz {
        //eprintln!("State {:x}", R[3 & THREE]);
        let mut c3: c_uchar =
            D[l0 as
              usize].R[Rans64DecGet(&mut R[3 & THREE],
                                    12 as uint32_t) as usize];
        out_buf[i4 as usize] = c3 as c_char;
        u8Rans64DecAdvanceSymbol(&mut R[THREE],
                                 pptr, &mut ptr_offset,
                                 &mut syms[l0 as usize][c3 as usize], 12i32 as uint32_t);
        //eprintln!("GET {:x} {:x}", l0, c3);
        //eprintln!("Atate {:x}", R[THREE]);
        l0 = c3 as c_int;
        i4 += 1
    }
    //*out_size = out_sz as c_uint;
}


#[no_mangle]
pub unsafe extern "C" fn rans_uncompress_O1b(mut in_0: *mut c_uchar,
                                             mut in_size: c_uint,
                                             mut out_size: *mut c_uint)
 -> *mut c_uchar {
    let mut cp: *mut c_uchar = in_0.offset(4isize);
    let mut ht_in_offset = 4usize;
    let mut i: c_int = 0;
    let mut j: c_int = -999i32;
    let mut x: c_int = 0;
    let mut out_sz: c_int = 0;
    let mut rle_i: c_int = 0;
    let mut rle_j: c_int = 0;
    let mut out_buf: *mut c_char = 0 as *mut c_char;
    let mut D: [ari_decoder; 256] = [ari_decoder{R: [0; 4096],}; 256];
    let mut syms: [[RansDecSymbol; 256]; 256] =
        [[RansDecSymbol{start: 0, freq: 0,}; 256]; 256];
    out_sz =
        (*in_0.offset(0isize) as c_int) << 0i32 |
            (*in_0.offset(1isize) as c_int) << 8i32 |
            (*in_0.offset(2isize) as c_int) << 16i32 |
            (*in_0.offset(3isize) as c_int) << 24i32;
    out_buf = malloc(out_sz as c_ulong) as *mut c_char;
    if out_buf.is_null() {
        return 0 as *mut c_uchar
    } else {
        rle_i = 0i32;
        let fresh70 = cp;
        cp = cp.offset(1);
        ht_in_offset += 1;
        i = *fresh70 as c_int;
        loop  {
            x = 0i32;
            rle_j = x;
            let fresh71 = cp;
            cp = cp.offset(1);
            ht_in_offset += 1;
            j = *fresh71 as c_int;
            loop  {
                let mut F: c_int = 0;
                let mut C: c_int = 0;
                let fresh72 = cp;
                cp = cp.offset(1);
                ht_in_offset += 1;
                F = *fresh72 as c_int;
                if F >= 128i32 {
                    F &= !128i32;
                    let fresh73 = cp;
                    cp = cp.offset(1);
                    ht_in_offset += 1;
                    F = (F & 127i32) << 8i32 | *fresh73 as c_int
                }
                C = x;
                if 0 == F { F = 1i32 << 12i32 }
                Rans64DecSymbolInit(&mut syms[i as usize][j as usize]
                                        , C as u16,
                                    F as u16);
                memset(&mut D[i as usize].R[x as usize] as *mut c_uchar
                           as *mut c_void, j, F as c_ulong);
                x += F;
                if 0 == rle_j && j + 1i32 == *cp as c_int {
                    let fresh74 = cp;
                    cp = cp.offset(1);
                    ht_in_offset += 1;
                    j = *fresh74 as c_int;
                    let fresh75 = cp;
                    cp = cp.offset(1);
                    ht_in_offset += 1;
                    rle_j = *fresh75 as c_int
                } else if 0 != rle_j {
                    rle_j -= 1;
                    j += 1
                } else {
                    let fresh76 = cp;
                    cp = cp.offset(1);
                    ht_in_offset += 1;
                    j = *fresh76 as c_int
                }
                if !(0 != j) { break ; }
            }
            if 0 == rle_i && i + 1i32 == *cp as c_int {
                let fresh77 = cp;
                cp = cp.offset(1);
                ht_in_offset += 1;
                i = *fresh77 as c_int;
                let fresh78 = cp;
                cp = cp.offset(1);
                ht_in_offset += 1;
                rle_i = *fresh78 as c_int
            } else if 0 != rle_i {
                rle_i -= 1;
                i += 1
            } else {
                let fresh79 = cp;
                cp = cp.offset(1);
                ht_in_offset += 1;
                i = *fresh79 as c_int
            }
            if !(0 != i) { break ; }
        }
        let mut rans0: RansState = 0;
        let mut rans1: RansState = 0;
        let mut rans2: RansState = 0;
        let mut rans3: RansState = 0;
        let mut ptr: *mut uint8_t = cp;
        let pptr = core::slice::from_raw_parts_mut(ptr as *mut u32, (in_size as usize - ht_in_offset) >> 2);
        let mut ptr_offset = 0usize;
        Rans64DecInit(&mut rans0,
                      pptr, &mut ptr_offset);
        Rans64DecInit(&mut rans1,
                      pptr, &mut ptr_offset);
        Rans64DecInit(&mut rans2,
                      pptr, &mut ptr_offset);
        Rans64DecInit(&mut rans3,
                      pptr, &mut ptr_offset);
        let mut isz4: c_int = out_sz >> 2i32;
        let mut l0: c_int = 0i32;
        let mut l1: c_int = 0i32;
        let mut l2: c_int = 0i32;
        let mut l3: c_int = 0i32;
        let mut i4: [c_int; 4] =
            [0i32 * isz4, 1i32 * isz4, 2i32 * isz4, 3i32 * isz4];
        let mut R: [RansState; 4] = [0; 4];
        R[0usize] = rans0;
        R[1usize] = rans1;
        R[2usize] = rans2;
        R[3usize] = rans3;
        while i4[0usize] < isz4 {
            let mut m: [uint32_t; 4] =
                [(R[0usize] &
                      (1u32 << 12i32).wrapping_sub(1i32 as c_uint) as
                          c_ulong) as uint32_t,
                 (R[1usize] &
                      (1u32 << 12i32).wrapping_sub(1i32 as c_uint) as
                          c_ulong) as uint32_t,
                 (R[2usize] &
                      (1u32 << 12i32).wrapping_sub(1i32 as c_uint) as
                          c_ulong) as uint32_t,
                 (R[3usize] &
                      (1u32 << 12i32).wrapping_sub(1i32 as c_uint) as
                          c_ulong) as uint32_t];
            let mut c: [uint8_t; 4] =
                [D[l0 as usize].R[m[0usize] as usize],
                 D[l1 as usize].R[m[1usize] as usize],
                 D[l2 as usize].R[m[2usize] as usize],
                 D[l3 as usize].R[m[3usize] as usize]];
            *out_buf.offset(i4[0usize] as isize) = c[0usize] as c_char;
            *out_buf.offset(i4[1usize] as isize) = c[1usize] as c_char;
            *out_buf.offset(i4[2usize] as isize) = c[2usize] as c_char;
            *out_buf.offset(i4[3usize] as isize) = c[3usize] as c_char;
            R[0usize] =
                (syms[l0 as usize][c[0usize] as usize].freq as
                     c_ulong).wrapping_mul(R[0usize] >> 12i32);
            R[0usize] =
                (R[0usize] as
                     c_ulong).wrapping_add(m[0usize].wrapping_sub(u32::from(syms[l0
                                                                                 as
                                                                                 usize][c[0usize]
                                                                                            as
                                                                                            usize].start))
                                                     as c_ulong) as
                    RansState as RansState;
            R[1usize] =
                (syms[l1 as usize][c[1usize] as usize].freq as
                     c_ulong).wrapping_mul(R[1usize] >> 12i32);
            R[1usize] =
                (R[1usize] as
                     c_ulong).wrapping_add(m[1usize].wrapping_sub(u32::from(syms[l1
                                                                                 as
                                                                                 usize][c[1usize]
                                                                                            as
                                                                                            usize].start))
                                                     as c_ulong) as
                    RansState as RansState;
            R[2usize] =
                (syms[l2 as usize][c[2usize] as usize].freq as
                     c_ulong).wrapping_mul(R[2usize] >> 12i32);
            R[2usize] =
                (R[2usize] as
                     c_ulong).wrapping_add(m[2usize].wrapping_sub(u32::from(syms[l2
                                                                                 as
                                                                                 usize][c[2usize]
                                                                                            as
                                                                                            usize].start))
                                                     as c_ulong) as
                    RansState as RansState;
            R[3usize] =
                (syms[l3 as usize][c[3usize] as usize].freq as
                     c_ulong).wrapping_mul(R[3usize] >> 12i32);
            R[3usize] =
                (R[3usize] as
                     c_ulong).wrapping_add(m[3usize].wrapping_sub(u32::from(syms[l3
                                                                                 as
                                                                                 usize][c[3usize]
                                                                                            as
                                                                                            usize].start))
                                                     as c_ulong) as
                    RansState as RansState;
            Rans64DecRenorm(&mut R[0usize],
                            pptr, &mut ptr_offset);
            Rans64DecRenorm(&mut R[1usize],
                            pptr, &mut ptr_offset);
            Rans64DecRenorm(&mut R[2usize],
                            pptr, &mut ptr_offset);
            Rans64DecRenorm(&mut R[3usize],
                            pptr, &mut ptr_offset);
            l0 = c[0usize] as c_int;
            l1 = c[1usize] as c_int;
            l2 = c[2usize] as c_int;
            l3 = c[3usize] as c_int;
            i4[0usize] += 1;
            i4[1usize] += 1;
            i4[2usize] += 1;
            i4[3usize] += 1
        }
        rans0 = R[0usize];
        rans1 = R[1usize];
        rans2 = R[2usize];
        rans3 = R[3usize];
        while i4[3usize] < out_sz {
            let mut c3: c_uchar =
                D[l3 as
                      usize].R[Rans64DecGet(&mut rans3,
                                            12i32 as uint32_t) as usize];
            *out_buf.offset(i4[3usize] as isize) = c3 as c_char;
            Rans64DecAdvanceSymbol(&mut rans3,
                                   pptr, &mut ptr_offset,
                                   &mut syms[l3 as usize][c3 as usize],
                                   12i32 as uint32_t);
            l3 = c3 as c_int;
            i4[3usize] += 1
        }
        *out_size = out_sz as c_uint;
        return out_buf as *mut c_uchar
    };
}
#[no_mangle]
pub unsafe extern "C" fn rans_compress(mut in_0: *mut c_uchar,
                                       mut in_size: c_uint,
                                       mut out_size: *mut c_uint,
                                       mut order: c_int)
 -> *mut c_uchar {
    return if 0 != order {
               rans_compress_O1d(in_0, in_size, out_size)
           } else { rans_compress_O0(in_0, in_size, out_size) };
}
#[no_mangle]
pub unsafe extern "C" fn rans_uncompress(mut in_0: *mut c_uchar,
                                         mut in_size: c_uint,
                                         mut out_size: *mut c_uint,
                                         mut order: c_int)
                                         -> *mut c_uchar {
    return if 0 != order {
               rans_uncompress_O1d(in_0, in_size, out_size)
           } else { rans_uncompress_O0b(in_0, in_size, out_size) };
}

