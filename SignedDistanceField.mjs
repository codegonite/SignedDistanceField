// https://shaderfun.com/2018/07/23/signed-distance-fields-part-8-gradients-bevels-and-noise/
// https://www.boristhebrave.com/2018/04/15/dual-contouring-tutorial/
// https://github.com/ssloy/least-squares-course
// https://www.mattkeeter.com/projects/contours/
// https://iquilezles.org/articles/distfunctions/

import fs from "node:fs";

export const EPSILON         = 1e-6;
export const INV_TWO_EPSILON = 1 / (2 * EPSILON);

function clamp(value, min, max) {
    return Math.min(max, Math.max(min, value));
}

function least_common_multiple(a, b) {
    return a * b / greatest_common_factor(a, b);
}

function greatest_common_factor(a, b) {
    if (Math.floor(a) != a || Math.floor(b) != b) {
        let numerator0   = a;
        let numerator1   = b;
        let denominator0 = 1;
        let denominator1 = 1;
        
        if (Math.floor(numerator0) != numerator0) {
            denominator0 = Math.pow(10, numerator0.toString().split('.')[1].length);
            numerator0 = Math.round(numerator0 * denominator0);
        }

        if (Math.floor(numerator1) != numerator1) {
            denominator1 = Math.pow(10, numerator1.toString().split('.')[1].length);
            numerator1 = Math.round(numerator1 * denominator1);
        }

        return greatest_common_factor(numerator0, numerator1) / least_common_multiple(denominator0, denominator1);
    }

    const remainder = a % b;
    if (remainder == 0) {
        return b;
    } else {
        return greatest_common_factor(b, remainder);
    }
}

export class Vector3 {
    static from(v) { return new Vector3(v.x, v.y, v.z); }
    static from_array(v) { return new Vector3(v[0], v[1], v[2]); }

    static neg(v0)                                  { return new Vector3(-v0.x,     -v0.y,      -v0.z); }
    static add(v0, v1)                              { return new Vector3(v0.x+v1.x,  v0.y+v1.y,  v0.z+v1.z); }
    static sub(v0, v1)                              { return new Vector3(v0.x-v1.x,  v0.y-v1.y,  v0.z-v1.z); }
    static mul(v0, v1)                              { return new Vector3(v0.x*v1.x,  v0.y*v1.y,  v0.z*v1.z); }
    static div(v0, v1)                              { return new Vector3(v0.x/v1.x,  v0.y/v1.y,  v0.z/v1.z); }
    static add_scalar(v0, s)                        { return new Vector3(v0.x+s,     v0.y+s,     v0.z+s); }
    static sub_scalar(v0, s)                        { return new Vector3(v0.x-s,     v0.y-s,     v0.z-s); }
    static mul_scalar(v0, s)                        { return new Vector3(v0.x*s,     v0.y*s,     v0.z*s); }
    static div_scalar(v0, s)                        { return new Vector3(v0.x/s,     v0.y/s,     v0.z/s); }
    static distance(v0, v1)                         { const x=v1.x-v0.x, y=v1.y-v0.y, z=v1.z-v0.z; return Math.sqrt(x*x + y*y + z*z); }
    static almost_equal(v0, v1)                     { const x=v1.x-v0.x, y=v1.y-v0.y, z=v1.z-v0.z; return Math.sqrt(x*x + y*y + z*z) <= EPSILON; }
    static magnitude_squared(v0)                    { const x=v0.x, y=v0.y, z=v0.z; return x*x + y*y + z*z; }
    static magnitude(v0)                            { const x=v0.x, y=v0.y, z=v0.z; return Math.sqrt(x*x + y*y + z*z); }
    static dot(v0, v1)                              { return v0.x*v1.x + v0.y*v1.y + v0.z*v1.z; }

    static average(vectors, out_vector = new Vector3()) {
        out_vector.set(0, 0, 0);
        for (let idx = 0; idx < vectors.length; ++idx) {
            out_vector.add(vectors[idx]);
        }
        out_vector.div_scalar(vectors.length);
        return out_vector;
    }

    static make_normal(a, b, c, out_vector = new Vector3()) {
        const ab0  = b.x - a.x;
        const ab1  = b.y - a.y;
        const ab2  = b.z - a.z;
        const ac0  = c.x - a.x;
        const ac1  = c.y - a.y;
        const ac2  = c.z - a.z;
        const out0 = ab1 * ac2 - ac1 * ab2;
        const out1 = ab2 * ac0 - ac2 * ab0;
        const out2 = ab0 * ac1 - ac0 * ab1;
        out_vector.set(out0, out1, out2);
        return out_vector.normalize();
    }

    static cross(v0, v1) {
        return new Vector3(v0.y * v1.z - v1.y * v0.z,
                           v0.z * v1.x - v1.z * v0.x,
                           v0.x * v1.y - v1.x * v0.y);
    }

    static normalize(v0) {
        const x0=v0.x, y0=v0.y, z0=v0.z;
        const length = Math.sqrt(x0*x0 + y0*y0 + z0*z0);

        if (length !== 0) {
            const factor = 1 / length;
            return new Vector3(x0 * factor, y0 * factor, z0 * factor);
        }

        return new Vector3(0, 0, 0);
    }

    static lerp(v0, v1, t0 = 0) {
        const x0 = v0.x, y0 = v0.y, z0 = v0.z;
        const x1 = v1.x, y1 = v1.y, z1 = v1.z;
        return new Vector3(t0 * (x1 - x0) + x0,
                           t0 * (y1 - y0) + y0,
                           t0 * (z1 - z0) + z0);
    }

    static transform(v, m) {
        const x = v.x, y = v.y, z = v.z;
        return new Vector3(x*m.x + y*m[4] + z*m[8]  + m[12],
                           x*m.y + y*m[5] + z*m[9]  + m[13],
                           x*m.z + y*m[6] + z*m[10] + m[14]);
    }

    constructor(x=0, y=0, z=0) {
        this.x = x;
        this.y = y;
        this.z = z;
    }

    get r() { return this.x; }   set r(v) { return this.x = v; }
    get g() { return this.y; }   set g(v) { return this.y = v; }
    get b() { return this.z; }   set b(v) { return this.z = v; }
    get [0]() { return this.x; } set [0](v) { return this.x = v; }
    get [1]() { return this.y; } set [1](v) { return this.y = v; }
    get [2]() { return this.z; } set [2](v) { return this.z = v; }
    *[Symbol.iterator]() { yield this.x; yield this.y; yield this.z; }

    /** Pure methods */
    clone()             { return new Vector3(this.x, this.y, this.z); }
    to_object()         { return { x: this.x, y: this.y, z: this.z }; }
    to_array()          { return [ this.x, this.y, this.z ]; }
    distance(other)     { const x=other.x-this.x, y=other.y-this.y, z=other.z-this.z; return Math.sqrt(x*x + y*y + z*z); }
    almost_equal(other) { const x=other.x-this.x, y=other.y-this.y, z=other.z-this.z; return Math.sqrt(x*x + y*y + z*z) <= EPSILON; }
    magnitude_squared() { return this.x*this.x + this.y*this.y + this.z*this.z; }
    magnitude()         { return Math.sqrt(this.x*this.x + this.y*this.y + this.z*this.z); }
    dot(other)          { return this.x*other.x + this.y*other.y + this.z*other.z; }
    lerp(other, t = 0) {
        const x0 = this.x,  y0 = this.y,  z0 = this.z;
        const x1 = other.x, y1 = other.y, z1 = other.z;

        return new Vector3((x1 - x0) * t + x0,
                           (y1 - y0) * t + y0,
                           (z1 - z0) * t + z0);
    }

    /** Modifying methods */
    abs()               { this.x = Math.abs(this.x); this.y = Math.abs(this.y); this.z = Math.abs(this.z); return this; }
    max(other)          { this.x = Math.max(this.x, other.x); this.y = Math.max(this.y, other.y); this.z = Math.max(this.z, other.z); return this; }
    min(other)          { this.x = Math.min(this.x, other.x); this.y = Math.min(this.y, other.y); this.z = Math.min(this.z, other.z); return this; }
    neg()               { this.x  = -this.x;  this.y = -this.y;  this.z = -this.z;  return this; }
    add(other)          { this.x +=  other.x; this.y += other.y; this.z += other.z; return this; }
    sub(other)          { this.x -=  other.x; this.y -= other.y; this.z -= other.z; return this; }
    mul(other)          { this.x *=  other.x; this.y *= other.y; this.z *= other.z; return this; }
    div(other)          { this.x /=  other.x; this.y /= other.y; this.z /= other.z; return this; }
    add_scalar(s)       { this.x +=  s;        this.y += s;        this.z += s;        return this; }
    sub_scalar(s)       { this.x -=  s;        this.y -= s;        this.z -= s;        return this; }
    mul_scalar(s)       { this.x *=  s;        this.y *= s;        this.z *= s;        return this; }
    div_scalar(s)       { this.x /=  s;        this.y /= s;        this.z /= s;        return this; }
    
    set_length(length) {
        const x0=this.x, y0=this.y, z0=this.z;
        const length_squared = x0*x0 + y0*y0 + z0*z0;

        if (length_squared == 0) {
            this.x = this.y = this.z = 0;
            return this;
        }

        const factor = length / Math.sqrt(length_squared);
        this.x = x0 * factor;
        this.y = y0 * factor;
        this.z = z0 * factor;
        return this;
    }

    normalize() {
        const x0=this.x, y0=this.y, z0=this.z;
        const length_squared = x0*x0 + y0*y0 + z0*z0;
        
        if (length_squared !== 0) {
            const factor = 1 / Math.sqrt(length_squared);
            this.x = x0 * factor;
            this.y = y0 * factor;
            this.z = z0 * factor;
            return this;
        }

        this.x = this.y = this.z = 0;
        return this;
    }

    transform(m) {
        const x = this.x, y = this.y, z = this.z;
        const m00=m.x,  m01=m.y,  m02=m.z;
        const m10=m[4],  m11=m[5],  m12=m[6];
        const m20=m[8],  m21=m[9],  m22=m[10];
        const m30=m[12], m31=m[13], m32=m[14];
        
        this.x = x*m00 + y*m10 + z*m20 + m30;
        this.y = x*m01 + y*m11 + z*m21 + m31;
        this.z = x*m02 + y*m12 + z*m22 + m32;
        return this;
    }
    
    cross(other) {
        const x0 = this.x,  y0 = this.y,  z0 = this.z;
        const x1 = other.x, y1 = other.y, z1 = other.z;
        this.x = y0 * z1 - y1 * z0;
        this.y = z0 * x1 - z1 * x0;
        this.z = x0 * y1 - x1 * y0;
        return this;
    }

    set(x, y, z)        { this.x =  x;          this.y =  y;          this.z =  z;         return this; }
    copy(other)         { this.x =  other.x;    this.y =  other.y;    this.z =  other.z;   return this; }
    neg_vector(other)   { this.x = -other.x;    this.y = -other.y;    this.z = -other.z;   return this; }
    add_vectors(v0, v1) { this.x =  v0.x+v1.x;  this.y =  v0.y+v1.y;  this.z =  v0.z+v1.z; return this; }
    sub_vectors(v0, v1) { this.x =  v0.x-v1.x;  this.y =  v0.y-v1.y;  this.z =  v0.z-v1.z; return this; }
    mul_vectors(v0, v1) { this.x =  v0.x*v1.x;  this.y =  v0.y*v1.y;  this.z =  v0.z*v1.z; return this; }
    div_vectors(v0, v1) { this.x =  v0.x/v1.x;  this.y =  v0.y/v1.y;  this.z =  v0.z/v1.z; return this; }

    lerp_vectors(v0, v1, t = 0) {
        const x0 = v0.x, y0 = v0.y, z0 = v0.z;
        const x1 = v1.x, y1 = v1.y, z1 = v1.z;

        this.x = (x1 - x0) * t + x0;
        this.y = (y1 - y0) * t + y0;
        this.z = (z1 - z0) * t + z0;
        return this;
    }

    transform_vector(m, v) {
        const x = v.x, y = v.y, z = v.z;
        const m00=m.x,  m01=m.y,  m02=m.z;
        const m10=m[4],  m11=m[5],  m12=m[6];
        const m20=m[8],  m21=m[9],  m22=m[10];
        const m30=m[12], m31=m[13], m32=m[14];
        
        this.x = x*m00 + y*m10 + z*m20 + m30;
        this.y = x*m01 + y*m11 + z*m21 + m31;
        this.z = x*m02 + y*m12 + z*m22 + m32;
        return this;
    }
    
    cross_vectors(v0, v1) {
        const x0 = v0.x, y0 = v0.y, z0 = v0.z;
        const x1 = v1.x, y1 = v1.y, z1 = v1.z;
        this.x = y0 * z1 - y1 * z0;
        this.y = z0 * x1 - z1 * x0;
        this.z = x0 * y1 - x1 * y0;
        return this;
    }
}

export class Matrix4 {
    static identity(out = new Matrix4()) {
        out.m00 = 1; out.m01 = 0; out.m02 = 0; out.m03 = 0;
        out.m10 = 0; out.m11 = 1; out.m12 = 0; out.m13 = 0;
        out.m20 = 0; out.m21 = 0; out.m22 = 1; out.m23 = 0;
        out.m30 = 0; out.m31 = 0; out.m32 = 0; out.m33 = 1;
        return out;
    }

    static scale(v) {
        const x=v[0], y=v[1], z=v[2];
        return new Matrix4(x, 0, 0, 0,
                           0, y, 0, 0,
                           0, 0, z, 0,
                           0, 0, 0, 1);
    }

    static translate(v) {
        const x=v[0], y=v[1], z=v[2];
        return new Matrix4(0, 0, 0, 0,
                           0, 0, 0, 0,
                           0, 0, 0, 0,
                           x, y, z, 1);
    }

    static rotatex(angle) {
        const c = Math.cos(angle);
        const s = Math.sin(angle);
        return new Matrix4(1, 0, 0, 0,
                           0, c,-s, 0,
                           0, s, c, 0,
                           0, 0, 0, 1);
    }

    static rotatey(angle) {
        const c = Math.cos(angle);
        const s = Math.sin(angle);
        return new Matrix4(c, 0, s, 0,
                           0, 1, 0, 0,
                          -s, 0, c, 0,
                           0, 0, 0, 1);
    }

    static rotatez(angle) {
        const c = Math.cos(angle);
        const s = Math.sin(angle);
        return new Matrix4(c,-s, 0, 0,
                           s, c, 0, 0,
                           0, 0, 1, 0,
                           0, 0, 0, 1);
    }

    static rotate(axis, angle) {
        const c = Math.cos(angle);
        const s = Math.sin(angle);
        const one_minus_c = 1 - c;

        let x=axis[0], y=axis[1], z=axis[2];
        const l = Math.sqrt(x * x + y * y + z * z);
        const f = 1 / l;
        
        x *= f;
        y *= f;
        z *= f;

        const xy = x*y, xs = x*s, xx = x*x;
        const xz = x*z, ys = y*s, yy = y*y;
        const yz = y*z, zs = z*s, zz = z*z;
        
        const b00 = c + xx*one_minus_c;
        const b10 = xy*one_minus_c - zs;
        const b20 = xz*one_minus_c + ys;
        
        const b01 = xy*one_minus_c + zs;
        const b11 = yy*one_minus_c + c;
        const b21 = yz*one_minus_c - xs;
        
        const b02 = xz*one_minus_c - ys;
        const b12 = yz*one_minus_c + xs;
        const b22 = zz*one_minus_c + c;

        return new Matrix4(b00, b10, b20, 0,
                           b01, b11, b21, 0,
                           b02, b12, b22, 0,
                           0,   0,   0,   1);
    }
    
    constructor(
        m00=1, m10=0, m20=0, m30=0,
        m01=0, m11=1, m21=0, m31=0,
        m02=0, m12=0, m22=1, m32=0,
        m03=0, m13=0, m23=0, m33=1,
    ) {
        this.m00 = m00;
        this.m01 = m01;
        this.m02 = m02;
        this.m03 = m03;
    
        this.m10 = m10;
        this.m11 = m11;
        this.m12 = m12;
        this.m13 = m13;
    
        this.m20 = m20;
        this.m21 = m21;
        this.m22 = m22;
        this.m23 = m23;
    
        this.m30 = m30;
        this.m31 = m31;
        this.m32 = m32;
        this.m33 = m33;
    }
    
    clone() {
        return new Matrix4(this.m00, this.m10, this.m20, this.m30,
                           this.m01, this.m11, this.m21, this.m31,
                           this.m02, this.m12, this.m22, this.m32,
                           this.m03, this.m13, this.m23, this.m33);
    }

    copy(other) {
        this.m00 = other.m00;
        this.m01 = other.m01;
        this.m02 = other.m02;
        this.m03 = other.m03;
        this.m10 = other.m10;
        this.m11 = other.m11;
        this.m12 = other.m12;
        this.m13 = other.m13;
        this.m20 = other.m20;
        this.m21 = other.m21;
        this.m22 = other.m22;
        this.m23 = other.m23;
        this.m30 = other.m30;
        this.m31 = other.m31;
        this.m32 = other.m32;
        this.m33 = other.m33;
        return this;
    }

    mul(other) {
        this.m00 = this.m00 * other.m00 + this.m10 * other.m01 + this.m20 * other.m02 + this.m30 * other.m03;
        this.m01 = this.m00 * other.m10 + this.m10 * other.m11 + this.m20 * other.m12 + this.m30 * other.m13;
        this.m02 = this.m00 * other.m20 + this.m10 * other.m21 + this.m20 * other.m22 + this.m30 * other.m23;
        this.m03 = this.m00 * other.m30 + this.m10 * other.m31 + this.m20 * other.m32 + this.m30 * other.m33;
        this.m10 = this.m01 * other.m00 + this.m11 * other.m01 + this.m21 * other.m02 + this.m31 * other.m03;
        this.m11 = this.m01 * other.m10 + this.m11 * other.m11 + this.m21 * other.m12 + this.m31 * other.m13;
        this.m12 = this.m01 * other.m20 + this.m11 * other.m21 + this.m21 * other.m22 + this.m31 * other.m23;
        this.m13 = this.m01 * other.m30 + this.m11 * other.m31 + this.m21 * other.m32 + this.m31 * other.m33;
        this.m20 = this.m02 * other.m00 + this.m12 * other.m01 + this.m22 * other.m02 + this.m32 * other.m03;
        this.m21 = this.m02 * other.m10 + this.m12 * other.m11 + this.m22 * other.m12 + this.m32 * other.m13;
        this.m22 = this.m02 * other.m20 + this.m12 * other.m21 + this.m22 * other.m22 + this.m32 * other.m23;
        this.m23 = this.m02 * other.m30 + this.m12 * other.m31 + this.m22 * other.m32 + this.m32 * other.m33;
        this.m30 = this.m03 * other.m00 + this.m13 * other.m01 + this.m23 * other.m02 + this.m33 * other.m03;
        this.m31 = this.m03 * other.m10 + this.m13 * other.m11 + this.m23 * other.m12 + this.m33 * other.m13;
        this.m32 = this.m03 * other.m20 + this.m13 * other.m21 + this.m23 * other.m22 + this.m33 * other.m23;
        this.m33 = this.m03 * other.m30 + this.m13 * other.m31 + this.m23 * other.m32 + this.m33 * other.m33;
        return this;
    }

    mul_matrices(m0, m1) {
        this.m00 = m0.m00 * m1.m00 + m0.m10 * m1.m01 + m0.m20 * m1.m02 + m0.m30 * m1.m03;
        this.m01 = m0.m00 * m1.m10 + m0.m10 * m1.m11 + m0.m20 * m1.m12 + m0.m30 * m1.m13;
        this.m02 = m0.m00 * m1.m20 + m0.m10 * m1.m21 + m0.m20 * m1.m22 + m0.m30 * m1.m23;
        this.m03 = m0.m00 * m1.m30 + m0.m10 * m1.m31 + m0.m20 * m1.m32 + m0.m30 * m1.m33;
        this.m10 = m0.m01 * m1.m00 + m0.m11 * m1.m01 + m0.m21 * m1.m02 + m0.m31 * m1.m03;
        this.m11 = m0.m01 * m1.m10 + m0.m11 * m1.m11 + m0.m21 * m1.m12 + m0.m31 * m1.m13;
        this.m12 = m0.m01 * m1.m20 + m0.m11 * m1.m21 + m0.m21 * m1.m22 + m0.m31 * m1.m23;
        this.m13 = m0.m01 * m1.m30 + m0.m11 * m1.m31 + m0.m21 * m1.m32 + m0.m31 * m1.m33;
        this.m20 = m0.m02 * m1.m00 + m0.m12 * m1.m01 + m0.m22 * m1.m02 + m0.m32 * m1.m03;
        this.m21 = m0.m02 * m1.m10 + m0.m12 * m1.m11 + m0.m22 * m1.m12 + m0.m32 * m1.m13;
        this.m22 = m0.m02 * m1.m20 + m0.m12 * m1.m21 + m0.m22 * m1.m22 + m0.m32 * m1.m23;
        this.m23 = m0.m02 * m1.m30 + m0.m12 * m1.m31 + m0.m22 * m1.m32 + m0.m32 * m1.m33;
        this.m30 = m0.m03 * m1.m00 + m0.m13 * m1.m01 + m0.m23 * m1.m02 + m0.m33 * m1.m03;
        this.m31 = m0.m03 * m1.m10 + m0.m13 * m1.m11 + m0.m23 * m1.m12 + m0.m33 * m1.m13;
        this.m32 = m0.m03 * m1.m20 + m0.m13 * m1.m21 + m0.m23 * m1.m22 + m0.m33 * m1.m23;
        this.m33 = m0.m03 * m1.m30 + m0.m13 * m1.m31 + m0.m23 * m1.m32 + m0.m33 * m1.m33;
        return this;
    }

    determinant() {
        const subfactor0 = this.m22 * this.m33 - this.m23 * this.m32;
        const subfactor1 = this.m21 * this.m33 - this.m23 * this.m31;
        const subfactor2 = this.m21 * this.m32 - this.m22 * this.m31;
        const subfactor3 = this.m20 * this.m33 - this.m23 * this.m30;
        const subfactor4 = this.m20 * this.m32 - this.m22 * this.m30;
        const subfactor5 = this.m20 * this.m31 - this.m21 * this.m30;

        const cofactor00 = this.m11 * subfactor0 - this.m12 * subfactor1 + this.m13 * subfactor2;
        const cofactor01 = this.m10 * subfactor0 - this.m12 * subfactor3 + this.m13 * subfactor4;
        const cofactor02 = this.m10 * subfactor1 - this.m11 * subfactor3 + this.m13 * subfactor5;
        const cofactor03 = this.m10 * subfactor2 - this.m11 * subfactor4 + this.m12 * subfactor5;

        return this.m00 * cofactor00 - this.m01 * cofactor01
             + this.m02 * cofactor02 - this.m03 * cofactor03;
    }

    transpose_inverse(m = this) {
        const subfactor0  = m.m22 * m.m33 - m.m23 * m.m32;
        const subfactor1  = m.m21 * m.m33 - m.m23 * m.m31;
        const subfactor2  = m.m21 * m.m32 - m.m22 * m.m31;
        const subfactor3  = m.m20 * m.m33 - m.m23 * m.m30;
        const subfactor4  = m.m20 * m.m32 - m.m22 * m.m30;
        const subfactor5  = m.m20 * m.m31 - m.m21 * m.m30;
        
        const subfactor6  = m.m02 * m.m13 - m.m03 * m.m12;
        const subfactor7  = m.m01 * m.m13 - m.m03 * m.m11;
        const subfactor8  = m.m01 * m.m12 - m.m02 * m.m11;
        const subfactor9  = m.m00 * m.m13 - m.m03 * m.m10;
        const subfactor10 = m.m00 * m.m12 - m.m02 * m.m10;
        const subfactor11 = m.m00 * m.m11 - m.m01 * m.m10;

        const cofactor00 = m.m11 * subfactor0  - m.m12 * subfactor1  + m.m13 * subfactor2;
        const cofactor01 = m.m12 * subfactor3  - m.m13 * subfactor4  - m.m10 * subfactor0;
        const cofactor02 = m.m10 * subfactor1  - m.m11 * subfactor3  + m.m13 * subfactor5;
        const cofactor03 = m.m11 * subfactor4  - m.m12 * subfactor5  - m.m10 * subfactor2;

        const cofactor10 = m.m02 * subfactor1  - m.m03 * subfactor2  - m.m01 * subfactor0;
        const cofactor11 = m.m00 * subfactor0  - m.m02 * subfactor3  + m.m03 * subfactor4;
        const cofactor12 = m.m01 * subfactor3  - m.m03 * subfactor5  - m.m00 * subfactor1;
        const cofactor13 = m.m00 * subfactor2  - m.m01 * subfactor4  + m.m02 * subfactor5;

        const cofactor20 = m.m31 * subfactor6  - m.m32 * subfactor7  + m.m33 * subfactor8;
        const cofactor21 = m.m32 * subfactor9  - m.m33 * subfactor10 - m.m30 * subfactor6;
        const cofactor22 = m.m30 * subfactor7  - m.m31 * subfactor9  + m.m33 * subfactor11;
        const cofactor23 = m.m31 * subfactor10 - m.m32 * subfactor11 - m.m30 * subfactor8;

        const cofactor30 = m.m22 * subfactor7  - m.m23 * subfactor8  - m.m21 * subfactor6;
        const cofactor31 = m.m20 * subfactor6  - m.m22 * subfactor9  + m.m23 * subfactor10;
        const cofactor32 = m.m21 * subfactor9  - m.m23 * subfactor11 - m.m20 * subfactor7;
        const cofactor33 = m.m20 * subfactor8  - m.m21 * subfactor10 + m.m22 * subfactor11;

        const determinant = m.m00 * cofactor00 + m.m01 * cofactor01 + m.m02 * cofactor02 + m.m03 * cofactor03;
        
        if (determinant === 0) {
            throw new Error("Determinant of the cofactor matrix cannot be zero!");
        }

        const factor = 1 / determinant;

        this.m00 = cofactor00*factor;
        this.m01 = cofactor01*factor;
        this.m02 = cofactor02*factor;
        this.m03 = cofactor03*factor;

        this.m10 = cofactor10*factor;
        this.m11 = cofactor11*factor;
        this.m12 = cofactor12*factor;
        this.m13 = cofactor13*factor;

        this.m20 = cofactor20*factor;
        this.m21 = cofactor21*factor;
        this.m22 = cofactor22*factor;
        this.m23 = cofactor23*factor;

        this.m30 = cofactor30*factor;
        this.m31 = cofactor31*factor;
        this.m32 = cofactor32*factor;
        this.m33 = cofactor33*factor;

        return this;
    }

    inverse(m = this) {
        const subfactor0  = m.m22 * m.m33 - m.m23 * m.m32;
        const subfactor1  = m.m21 * m.m33 - m.m23 * m.m31;
        const subfactor2  = m.m21 * m.m32 - m.m22 * m.m31;
        const subfactor3  = m.m20 * m.m33 - m.m23 * m.m30;
        const subfactor4  = m.m20 * m.m32 - m.m22 * m.m30;
        const subfactor5  = m.m20 * m.m31 - m.m21 * m.m30;
        
        const subfactor6  = m.m02 * m.m13 - m.m03 * m.m12;
        const subfactor7  = m.m01 * m.m13 - m.m03 * m.m11;
        const subfactor8  = m.m01 * m.m12 - m.m02 * m.m11;
        const subfactor9  = m.m00 * m.m13 - m.m03 * m.m10;
        const subfactor10 = m.m00 * m.m12 - m.m02 * m.m10;
        const subfactor11 = m.m00 * m.m11 - m.m01 * m.m10;

        const cofactor00 = m.m11 * subfactor0  - m.m12 * subfactor1  + m.m13 * subfactor2;
        const cofactor01 = m.m12 * subfactor3  - m.m13 * subfactor4  - m.m10 * subfactor0;
        const cofactor02 = m.m10 * subfactor1  - m.m11 * subfactor3  + m.m13 * subfactor5;
        const cofactor03 = m.m11 * subfactor4  - m.m12 * subfactor5  - m.m10 * subfactor2;

        const cofactor10 = m.m02 * subfactor1  - m.m03 * subfactor2  - m.m01 * subfactor0;
        const cofactor11 = m.m00 * subfactor0  - m.m02 * subfactor3  + m.m03 * subfactor4;
        const cofactor12 = m.m01 * subfactor3  - m.m03 * subfactor5  - m.m00 * subfactor1;
        const cofactor13 = m.m00 * subfactor2  - m.m01 * subfactor4  + m.m02 * subfactor5;

        const cofactor20 = m.m31 * subfactor6  - m.m32 * subfactor7  + m.m33 * subfactor8;
        const cofactor21 = m.m32 * subfactor9  - m.m33 * subfactor10 - m.m30 * subfactor6;
        const cofactor22 = m.m30 * subfactor7  - m.m31 * subfactor9  + m.m33 * subfactor11;
        const cofactor23 = m.m31 * subfactor10 - m.m32 * subfactor11 - m.m30 * subfactor8;

        const cofactor30 = m.m22 * subfactor7  - m.m23 * subfactor8  - m.m21 * subfactor6;
        const cofactor31 = m.m20 * subfactor6  - m.m22 * subfactor9  + m.m23 * subfactor10;
        const cofactor32 = m.m21 * subfactor9  - m.m23 * subfactor11 - m.m20 * subfactor7;
        const cofactor33 = m.m20 * subfactor8  - m.m21 * subfactor10 + m.m22 * subfactor11;

        const determinant = m.m00 * cofactor00 + m.m01 * cofactor01 + m.m02 * cofactor02 + m.m03 * cofactor03;
        
        if (determinant === 0) {
            throw new Error("Determinant of the cofactor matrix cannot be zero!");
        }

        const factor = 1 / determinant;

        this.m00 = cofactor00*factor;
        this.m01 = cofactor10*factor;
        this.m02 = cofactor20*factor;
        this.m03 = cofactor30*factor;

        this.m10 = cofactor01*factor;
        this.m11 = cofactor11*factor;
        this.m12 = cofactor21*factor;
        this.m13 = cofactor31*factor;

        this.m20 = cofactor02*factor;
        this.m21 = cofactor12*factor;
        this.m22 = cofactor22*factor;
        this.m23 = cofactor32*factor;

        this.m30 = cofactor03*factor;
        this.m31 = cofactor13*factor;
        this.m32 = cofactor23*factor;
        this.m33 = cofactor33*factor;

        return this;
    }

    is_identity() {
        return this.m00 === 1 && this.m01 === 0 && this.m02 === 0 && this.m03 === 0
        &&     this.m10 === 0 && this.m11 === 1 && this.m12 === 0 && this.m13 === 0
        &&     this.m20 === 0 && this.m21 === 0 && this.m22 === 1 && this.m23 === 0
        &&     this.m30 === 0 && this.m31 === 0 && this.m32 === 0 && this.m33 === 1;
    }

    equals(other) {
        return Math.abs(this.m00 - other.m00) <= EPSILON && Math.abs(this.m10 - other.m10) <= EPSILON
        &&     Math.abs(this.m20 - other.m20) <= EPSILON && Math.abs(this.m30 - other.m30) <= EPSILON
        &&     Math.abs(this.m01 - other.m01) <= EPSILON && Math.abs(this.m11 - other.m11) <= EPSILON
        &&     Math.abs(this.m21 - other.m21) <= EPSILON && Math.abs(this.m31 - other.m31) <= EPSILON
        &&     Math.abs(this.m02 - other.m02) <= EPSILON && Math.abs(this.m12 - other.m12) <= EPSILON
        &&     Math.abs(this.m22 - other.m22) <= EPSILON && Math.abs(this.m32 - other.m32) <= EPSILON
        &&     Math.abs(this.m03 - other.m03) <= EPSILON && Math.abs(this.m13 - other.m13) <= EPSILON
        &&     Math.abs(this.m23 - other.m23) <= EPSILON && Math.abs(this.m33 - other.m33) <= EPSILON;
    }

    scale(v, m=this) {
        this.m00 = m.m00*v.x;
        this.m01 = m.m01*v.y;
        this.m02 = m.m02*v.z;
        this.m10 = m.m10*v.x;
        this.m11 = m.m11*v.y;
        this.m12 = m.m12*v.z;
        this.m20 = m.m20*v.x;
        this.m21 = m.m21*v.y;
        this.m22 = m.m22*v.z;
        this.m30 = m.m30*v.x;
        this.m31 = m.m31*v.y;
        this.m32 = m.m32*v.z;
        return this;
    }

    translate(v, m=this) {
        this.m03 = m.m00*v.x + m.m01*v.y + m.m02*v.z + m.m03;
        this.m13 = m.m10*v.x + m.m11*v.y + m.m12*v.z + m.m13;
        this.m23 = m.m20*v.x + m.m21*v.y + m.m22*v.z + m.m23;
        this.m33 = m.m30*v.x + m.m31*v.y + m.m32*v.z + m.m33;
        return this;
    }

    rotatex(angle, m=this) {
        const c = Math.cos(angle);
        const s = Math.sin(angle);
        
        this.m01 = m.m01*c - m.m02*s;
        this.m02 = m.m02*c + m.m01*s;
        this.m11 = m.m11*c - m.m12*s;
        this.m12 = m.m12*c + m.m11*s;
        this.m21 = m.m21*c - m.m22*s;
        this.m22 = m.m22*c + m.m21*s;
        this.m31 = m.m31*c - m.m32*s;
        this.m32 = m.m32*c + m.m31*s;
        return this;
    }

    rotatey(angle, m=this) {
        const c = Math.cos(angle);
        const s = Math.sin(angle);
        
        this.m00 = m.m00*c + m.m02*s;
        this.m02 = m.m02*c - m.m00*s;
        this.m10 = m.m10*c + m.m12*s;
        this.m12 = m.m12*c - m.m10*s;
        this.m20 = m.m20*c + m.m22*s;
        this.m22 = m.m22*c - m.m20*s;
        this.m30 = m.m30*c + m.m32*s;
        this.m32 = m.m32*c - m.m30*s;
        return this;
    }

    rotatez(angle, m=this) {
        const c = Math.cos(angle);
        const s = Math.sin(angle);

        this.m00 = m.m00*c - m.m01*s;
        this.m01 = m.m01*c + m.m00*s;
        this.m10 = m.m10*c - m.m11*s;
        this.m11 = m.m11*c + m.m10*s;
        this.m20 = m.m20*c - m.m21*s;
        this.m21 = m.m21*c + m.m20*s;
        this.m30 = m.m30*c - m.m31*s;
        this.m31 = m.m31*c + m.m30*s;
        return this;
    }

    rotate(axis, angle, m=this) {
        const c = Math.cos(angle);
        const s = Math.sin(angle);
        const one_minus_c = 1 - c;

        let x=axis[0], y=axis[1], z=axis[2];
        const l = Math.sqrt(x * x + y * y + z * z);
        const f = 1 / l;
        
        x *= f;
        y *= f;
        z *= f;

        const xy = x*y, xs = x*s, xx = x*x;
        const xz = x*z, ys = y*s, yy = y*y;
        const yz = y*z, zs = z*s, zz = z*z;
        
        const b00 = c + xx*one_minus_c;
        const b10 = xy*one_minus_c - zs;
        const b20 = xz*one_minus_c + ys;
        
        const b01 = xy*one_minus_c + zs;
        const b11 = yy*one_minus_c + c;
        const b21 = yz*one_minus_c - xs;
        
        const b02 = xz*one_minus_c - ys;
        const b12 = yz*one_minus_c + xs;
        const b22 = zz*one_minus_c + c;

        this.m00 = m.m10 * b10 + m.m20 * b20 + m.m00 * b00; // m00
        this.m01 = m.m00 * b01 + m.m20 * b21 + m.m10 * b11; // m01
        this.m02 = m.m00 * b02 + m.m10 * b12 + m.m20 * b22; // m02
        this.m10 = m.m11 * b10 + m.m21 * b20 + m.m01 * b00; // m10
        this.m11 = m.m01 * b01 + m.m21 * b21 + m.m11 * b11; // m11
        this.m12 = m.m01 * b02 + m.m11 * b12 + m.m21 * b22; // m12
        this.m20 = m.m12 * b10 + m.m22 * b20 + m.m02 * b00; // m20
        this.m21 = m.m02 * b01 + m.m22 * b21 + m.m12 * b11; // m21
        this.m22 = m.m02 * b02 + m.m12 * b12 + m.m22 * b22; // m22
        this.m30 = m.m13 * b10 + m.m23 * b20 + m.m03 * b00; // m30
        this.m31 = m.m03 * b01 + m.m23 * b21 + m.m13 * b11; // m31
        this.m32 = m.m03 * b02 + m.m13 * b12 + m.m23 * b22; // m32
        return this;
    }
}

export class Polygon {
    constructor(positions) {
        this.positions = positions;
    }
}

export class Mesh {
    constructor(polygons) {
        this.polygons = polygons;
    }

    serialize_to_obj(filepath) {
        const stream = fs.createWriteStream(filepath);

        for (let idx = 0, position_idx = 1; idx < this.polygons.length; ++idx) {
            const polygon = this.polygons[idx];
            for (const vertex of polygon.positions) {
                stream.write(`v\t${vertex.x}\t${vertex.y}\t${vertex.z}\n`);
            }

            if (polygon.positions.length != 0) {
                stream.write(`f\t${polygon.positions.map(()=>position_idx++).join("\t")}\n`);
            }
        }

        stream.close(err => {
            console.log(`Written ${stream.bytesWritten} bytes to file "${filepath}" in Wavefront OBJ file format!`);
        });
    }
}

export class SamplePoint {
    constructor(x, y, z, distance) {
        this.x = x;
        this.y = y;
        this.z = z;
        this.distance = distance;
    }

    get_intersection(other, field, out = new Vector3()) {
        if (this.distance * other.distance > 0) {
            return null;
        }

        const PRECISION = 0.00000005;
        const dx = other.x - this.x;
        const dy = other.y - this.y;
        const dz = other.z - this.z;
        const total_distance = Math.sqrt(dx*dx + dy*dy + dz*dz);
        const inv_total_distance = 1 / total_distance;
        const norm_dx = dx * inv_total_distance;
        const norm_dy = dy * inv_total_distance;
        const norm_dz = dz * inv_total_distance;
        let sphere_distance = 0;

        out.set(this.x, this.y, this.z);
        do {
            sphere_distance = field.calculate_signed_distance(out.x, out.y, out.z);
            out.x += norm_dx * sphere_distance;
            out.y += norm_dy * sphere_distance;
            out.z += norm_dz * sphere_distance;
        } while (sphere_distance > PRECISION);

        return out;
    }

    interpolate_distance(other, out = new Vector3()) {
        if (this.distance * other.distance > 0) {
            return null;
        }

        const total_distance = Math.abs(this.distance) + Math.abs(other.distance);

        if (total_distance == 0) {
            return null;
        }

        return out.lerp_vectors(this, other, Math.abs(this.distance) / total_distance);
    }
};

export class SignedDistanceField {
    static temp_vector0 = new Vector3();

    bounding_box = new BoundingBox();
    
    calculate_signed_distance(x, y, z) {
        throw new Error("calculate_signed_distance is unimplemented!");
    }

    calculate_gradient(x, y, z, out_vector = new Vector3(0, 0, 0)) {
        const distance_max_x = this.calculate_signed_distance(x + EPSILON, y, z);
        const distance_min_x = this.calculate_signed_distance(x - EPSILON, y, z);
        const distance_max_y = this.calculate_signed_distance(x, y + EPSILON, z);
        const distance_min_y = this.calculate_signed_distance(x, y - EPSILON, z);
        const distance_max_z = this.calculate_signed_distance(x, y, z + EPSILON);
        const distance_min_z = this.calculate_signed_distance(x, y, z - EPSILON);
        out_vector.set((distance_max_x - distance_min_x) * INV_TWO_EPSILON,
                       (distance_max_y - distance_min_y) * INV_TWO_EPSILON,
                       (distance_max_z - distance_min_z) * INV_TWO_EPSILON);
        return out_vector;
    }

    calculate_to_surface_vector(vector, out_vector = new Vector3(0, 0, 0)) {
        const x = vector.x, y = vector.y, z = vector.z;
        const distance = this.calculate_signed_distance(x, y, z);
        return this.calculate_gradient(x, y, z, out_vector).set_length(distance);
    }

    calculate_closest_surface_vector(x, y, z, out_vector = new Vector3(0, 0, 0)) {
        const distance = this.calculate_signed_distance(x, y, z);
        const delta = this.calculate_gradient(x, y, z, out_vector).set_length(distance);
        console.log({x, y, z, distance});
        return delta.add(SignedDistanceField.temp_vector0.set(x, y, z));
    }

    get_sample_point(x, y, z) {
        const d = this.calculate_signed_distance(x, y, z);
        return new SamplePoint(x, y, z, d);
    }

    calculate_bounding_box() {
        const LARGEST_NUMBER = 1000000;
        const min_x = this.calculate_closest_surface_vector(-LARGEST_NUMBER, 0, 0).x;
        const max_x = this.calculate_closest_surface_vector( LARGEST_NUMBER, 0, 0).x;
        const min_y = this.calculate_closest_surface_vector(0, -LARGEST_NUMBER, 0).y;
        const max_y = this.calculate_closest_surface_vector(0,  LARGEST_NUMBER, 0).y;
        const min_z = this.calculate_closest_surface_vector(0, 0, -LARGEST_NUMBER).z;
        const max_z = this.calculate_closest_surface_vector(0, 0,  LARGEST_NUMBER).z;
        return new BoundingBox(min_x, max_x, min_y, max_y, min_z, max_z);
    }

    union(... others)                                  { return new SignedDistanceFieldUnion(this, ... others); }
    difference(... others)                             { return new SignedDistanceFieldDifference(this, ... others); }
    intersection(... others)                           { return new SignedDistanceFieldIntersection(this, ... others); }
    union_smooth(smoothness, ... others)               { return new SignedDistanceFieldUnionSmooth(smoothness, this, ... others); }
    difference_smooth(smoothness, ... others)          { return new SignedDistanceFieldDifferenceSmooth(smoothness, this, ... others); }
    intersection_smooth(smoothness, ... others)        { return new SignedDistanceFieldIntersectionSmooth(smoothness, this, ... others); }
    transform(m)                                       { return new SignedDistanceFieldTransform(this, m); }
    translate(x, y, z)                                 { return new SignedDistanceFieldTranslate(x, y, z, this); }
    offset(radius)                                     { return new SignedDistanceFieldOffset(this, radius); }

    surface_nets(step_count = 100) {
        const faces = [];
        const sample_points = [];
        const step_size = this.bounding_box.largest_side / step_count;
        console.log(this.bounding_box.min_x);
        console.log(this.bounding_box.min_y);
        console.log(this.bounding_box.min_z);
        console.log(this.bounding_box.max_x);
        console.log(this.bounding_box.max_y);
        console.log(this.bounding_box.max_z);
        console.log(this.bounding_box.size_x);
        console.log(this.bounding_box.size_y);
        console.log(this.bounding_box.size_z);
        console.log(this.bounding_box.largest_side)
        console.log({step_count});
        console.log({step_size});
        
        const smallest_sample_count = this.bounding_box.smallest_side / this.bounding_box.smallest_cube_size();
        const bounding_box          = smallest_sample_count >= 1000 ? this.bounding_box.clone().optimize_smallest_cube_size(step_size) : this.bounding_box;
        const smallest_cube_size    = bounding_box.smallest_cube_size();

        const cube_size = step_size >= smallest_cube_size ? smallest_cube_size / Math.floor(Math.max(smallest_cube_size / step_size, 1))
                                                          : step_size / Math.floor(Math.max(step_size / smallest_cube_size, 1));
        const sample_count_x = Math.round(bounding_box.size_x / cube_size) + 2;
        const sample_count_y = Math.round(bounding_box.size_y / cube_size) + 2;
        const sample_count_z = Math.round(bounding_box.size_z / cube_size) + 2;
        const half_size      = cube_size / 2;
        const axis0_stride   = sample_count_y * sample_count_z;
        const axis1_stride   = sample_count_z;

        for (let idx0 = 0, x = bounding_box.min_x - half_size; idx0 < sample_count_x; ++idx0, x += cube_size) {
            for (let idx1 = 0, y = bounding_box.min_y - half_size; idx1 < sample_count_y; ++idx1, y += cube_size) {
                for (let idx2 = 0, z = bounding_box.min_z - half_size; idx2 < sample_count_z; ++idx2, z += cube_size) {
                    const sample_point = this.get_sample_point(x, y, z);
                    sample_points.push(sample_point);
                }
            }
        }

        let index = 0;
        for (let idx0 = 0, x = bounding_box.min_x - half_size; idx0 < sample_count_x; ++idx0, x += cube_size) {
            for (let idx1 = 0, y = bounding_box.min_y - half_size; idx1 < sample_count_y; ++idx1, y += cube_size) {
                for (let idx2 = 0, z = bounding_box.min_z - half_size; idx2 < sample_count_z; ++idx2, z += cube_size, ++index) {
                    if (sample_points[index].distance > 0) {
                        continue;
                    }

                    const top_index    = index + 1;
                    const bottom_index = index - 1;
                    const back_index   = index - axis1_stride;
                    const front_index  = index + axis1_stride;
                    const right_index  = index + axis0_stride;
                    const left_index   = index - axis0_stride;

                    if (top_index >= 0 && top_index < sample_points.length
                    &&  sample_points[index].distance * sample_points[top_index].distance <= 0) {
                        faces.push(new Polygon([
                            new Vector3(x - half_size, y - half_size, z + half_size),
                            new Vector3(x + half_size, y - half_size, z + half_size),
                            new Vector3(x + half_size, y + half_size, z + half_size),
                            new Vector3(x - half_size, y + half_size, z + half_size),
                        ]));
                    }

                    if (bottom_index >= 0 && bottom_index < sample_points.length
                    &&  sample_points[index].distance * sample_points[bottom_index].distance <= 0) {
                        faces.push(new Polygon([
                            new Vector3(x - half_size, y + half_size, z - half_size),
                            new Vector3(x + half_size, y + half_size, z - half_size),
                            new Vector3(x + half_size, y - half_size, z - half_size),
                            new Vector3(x - half_size, y - half_size, z - half_size),
                        ]));
                    }

                    if (back_index >= 0 && back_index < sample_points.length
                    &&  sample_points[index].distance * sample_points[back_index].distance <= 0) {
                        faces.push(new Polygon([
                            new Vector3(x - half_size, y - half_size, z - half_size),
                            new Vector3(x + half_size, y - half_size, z - half_size),
                            new Vector3(x + half_size, y - half_size, z + half_size),
                            new Vector3(x - half_size, y - half_size, z + half_size),
                        ]));
                    }

                    if (front_index >= 0 && front_index < sample_points.length
                    &&  sample_points[index].distance * sample_points[front_index].distance <= 0) {
                        faces.push(new Polygon([
                            new Vector3(x - half_size, y + half_size, z + half_size),
                            new Vector3(x + half_size, y + half_size, z + half_size),
                            new Vector3(x + half_size, y + half_size, z - half_size),
                            new Vector3(x - half_size, y + half_size, z - half_size),
                        ]));
                    }

                    if (right_index >= 0 && right_index < sample_points.length
                    &&  sample_points[index].distance * sample_points[right_index].distance <= 0) {
                        faces.push(new Polygon([
                            new Vector3(x + half_size, y - half_size, z - half_size),
                            new Vector3(x + half_size, y + half_size, z - half_size),
                            new Vector3(x + half_size, y + half_size, z + half_size),
                            new Vector3(x + half_size, y - half_size, z + half_size),
                        ]));
                    }

                    if (left_index >= 0 && left_index < sample_points.length
                    &&  sample_points[index].distance * sample_points[left_index].distance <= 0) {
                        faces.push(new Polygon([
                            new Vector3(x - half_size, y - half_size, z + half_size),
                            new Vector3(x - half_size, y + half_size, z + half_size),
                            new Vector3(x - half_size, y + half_size, z - half_size),
                            new Vector3(x - half_size, y - half_size, z - half_size),
                        ]));
                    }
                }
            }
        }

        for (let idx0 = 0; idx0 < faces.length; ++idx0) {
            const face = faces[idx0];
            for (let idx1 = 0; idx1 < faces[idx1].positions.length; ++idx1) {
                face.positions[idx1].sub(this.calculate_to_surface_vector(face.positions[idx1]));
                face.positions[idx1].sub(this.calculate_to_surface_vector(face.positions[idx1]));
                face.positions[idx1].sub(this.calculate_to_surface_vector(face.positions[idx1]));
                face.positions[idx1].sub(this.calculate_to_surface_vector(face.positions[idx1]));
            }
        }

        return new Mesh(faces);
    }
}

export class BoundingBox {
    get size_x() { return this.max_x - this.min_x; }
    get size_y() { return this.max_y - this.min_y; }
    get size_z() { return this.max_z - this.min_z; }
    get center_x() { return (this.max_x + this.min_x) / 2; }
    get center_y() { return (this.max_y + this.min_y) / 2; }
    get center_z() { return (this.max_z + this.min_z) / 2; }
    get largest_side() { return Math.max(this.size_x, this.size_y, this.size_z); }
    get smallest_side() { return Math.min(this.size_x, this.size_y, this.size_z); }

    constructor(min_x = 0, min_y = 0, min_z = 0, max_x = 0, max_y = 0, max_z = 0) {
        this.min_x = min_x;
        this.min_y = min_y;
        this.min_z = min_z;
        this.max_x = max_x;
        this.max_y = max_y;
        this.max_z = max_z;
    }

    static enclose_distance_fields(fields) {
        const bounding_box = fields[0].bounding_box.clone();
        
        for (let idx = 1; idx < fields.length; ++idx) {
            bounding_box.enclose_box(fields[idx].bounding_box);
        }

        return bounding_box;
    }

    static intersecting_distance_fields(fields) {
        const bounding_box = fields[0].bounding_box.clone();
        
        for (let idx = 1; idx < fields.length; ++idx) {
            bounding_box.intersect_box(fields[idx].bounding_box);
        }

        return bounding_box;
    }

    static create_enclosing_box(box0, box1) {
        return new BoundingBox(
            Math.min(box0.min_x, box1.min_x),
            Math.min(box0.min_y, box1.min_y),
            Math.min(box0.min_z, box1.min_z),
            Math.max(box0.max_x, box1.max_x),
            Math.max(box0.max_y, box1.max_y),
            Math.max(box0.max_z, box1.max_z),
        );
    }

    static create_intersecting_box(box0, box1) {
        return new BoundingBox(
            Math.max(box0.min_x, box1.min_x),
            Math.max(box0.min_y, box1.min_y),
            Math.max(box0.min_z, box1.min_z),
            Math.min(box0.max_x, box1.max_x),
            Math.min(box0.max_y, box1.max_y),
            Math.min(box0.max_z, box1.max_z),
        );
    }
    
    center() {
        return new Vector3(this.center_x, this.center_y, this.center_z);
    }

    smallest_cube_size() {
        return greatest_common_factor(this.size_x,
               greatest_common_factor(this.size_y,
                                      this.size_z));
    }

    optimize_smallest_cube_size(step_size) {
        const center_x = (this.max_x + this.min_x) / 2;
        const center_y = (this.max_y + this.min_y) / 2;
        const center_z = (this.max_z + this.min_z) / 2;
        const size_x = step_size * Math.round(this.size_x / step_size);
        const size_y = step_size * Math.round(this.size_y / step_size);
        const size_z = step_size * Math.round(this.size_z / step_size);

        this.min_x = center_x - size_x / 2;
        this.min_y = center_y - size_y / 2;
        this.min_z = center_z - size_z / 2;
        this.max_x = center_x + size_x / 2;
        this.max_y = center_y + size_y / 2;
        this.max_z = center_z + size_z / 2;

        return this;
    }

    clone() {
        return new BoundingBox(
            this.min_x,
            this.min_y,
            this.min_z,
            this.max_x,
            this.max_y,
            this.max_z,
        );
    }

    grow(size) {
        this.min_x -= size;
        this.min_y -= size;
        this.min_z -= size;
        this.max_x += size;
        this.max_y += size;
        this.max_z += size;
        return this;
    }

    mul(factor) {
        this.min_x *= factor;
        this.min_y *= factor;
        this.min_z *= factor;
        this.max_x *= factor;
        this.max_y *= factor;
        this.max_z *= factor;
        return this;
    }

    to_polygons(polygons = []) {
        const position0 = new Vector3(this.min_x, this.min_y, this.min_z);
        const position1 = new Vector3(this.max_x, this.min_y, this.min_z);
        const position2 = new Vector3(this.max_x, this.max_y, this.min_z);
        const position3 = new Vector3(this.min_x, this.max_y, this.min_z);
        const position4 = new Vector3(this.min_x, this.min_y, this.max_z);
        const position5 = new Vector3(this.max_x, this.min_y, this.max_z);
        const position6 = new Vector3(this.max_x, this.max_y, this.max_z);
        const position7 = new Vector3(this.min_x, this.max_y, this.max_z);
        polygons.push(new Polygon([ position3, position2, position1, position0 ]));
        polygons.push(new Polygon([ position4, position5, position6, position7 ]));
        polygons.push(new Polygon([ position0, position1, position5, position4 ]));
        polygons.push(new Polygon([ position2, position3, position7, position6 ]));
        polygons.push(new Polygon([ position1, position2, position6, position5 ]));
        polygons.push(new Polygon([ position4, position7, position3, position0 ]));
        return polygons;
    }

    enclose_box(other) {
        this.min_x = Math.min(this.min_x, other.min_x);
        this.min_y = Math.min(this.min_y, other.min_y);
        this.min_z = Math.min(this.min_z, other.min_z);
        this.max_x = Math.max(this.max_x, other.max_x);
        this.max_y = Math.max(this.max_y, other.max_y);
        this.max_z = Math.max(this.max_z, other.max_z);
        return this;
    }

    intersect_box(other) {
        this.min_x = Math.max(this.min_x, other.min_x);
        this.min_y = Math.max(this.min_y, other.min_y);
        this.min_z = Math.max(this.min_z, other.min_z);
        this.max_x = Math.min(this.max_x, other.max_x);
        this.max_y = Math.min(this.max_y, other.max_y);
        this.max_z = Math.min(this.max_z, other.max_z);
        return this;
    }

}

export class SignedDistanceFieldSphere extends SignedDistanceField {
    constructor(position = new Vector3(0, 0, 0), radius = 5) {
        super();
        this.position_x = position[0];
        this.position_y = position[1];
        this.position_z = position[2];
        this.radius     = radius;
        this.bounding_box = new BoundingBox(
            this.position_x - radius,
            this.position_y - radius,
            this.position_z - radius,
            this.position_x + radius,
            this.position_y + radius,
            this.position_z + radius,
        );
    }

    calculate_signed_distance(x, y, z) {
        const dx = x - this.position_x;
        const dy = y - this.position_y;
        const dz = z - this.position_z;
        return Math.sqrt(dx*dx + dy*dy + dz*dz) - this.radius;
    }
}

export class SignedDistanceFieldCylinder extends SignedDistanceField {
    constructor(position = new Vector3(0, 0, 0), height = 20, radius = 10) {
        super();
        this.position_x  = position[0];
        this.position_y  = position[1];
        this.position_z  = position[2];
        this.radius      = radius;
        this.half_height = height / 2;
        this.bounding_box = new BoundingBox(
            this.position_x - this.radius,
            this.position_y - this.radius,
            this.position_z - this.half_height,
            this.position_x + this.radius,
            this.position_y + this.radius,
            this.position_z + this.half_height,
        );
    }

    calculate_signed_distance(x, y, z) {
        const dx = x - this.position_x;
        const dy = y - this.position_y;
        const dz = z - this.position_z;
        const d0 = Math.abs(Math.sqrt(dx*dx + dy*dy)) - this.radius;
        const d1 = Math.abs(dz) - this.half_height;
        const a0 = Math.max(d0,0.0);
        const a1 = Math.max(d1,0.0);
        return Math.min(Math.max(d0,d1),0.0) + Math.sqrt(a0*a0 + a1*a1);
    }
}

export class SignedDistanceFieldCube extends SignedDistanceField {
    constructor(position = [0, 0, 0], size = 10, corner_radius = 0) {
        if (typeof size != "number") {
            throw new Error("Invalid 'size' argument must be a number!");
        }

        if (typeof corner_radius != "number") {
            throw new Error("Invalid 'corner_radius' argument must be a number!");
        }

        if (position == undefined || typeof position[0] != "number" || typeof position[1] != "number" || typeof position[2] != "number") {
            throw new Error("Invalid 'position' argument must be a number array!");
        }

        super();
        this.position_x    = position[0];
        this.position_y    = position[1];
        this.position_z    = position[2];
        this.half_size     = size / 2;
        this.corner_radius = corner_radius;
        this.bounding_box = new BoundingBox(
            this.position_x - this.half_size,
            this.position_y - this.half_size,
            this.position_z - this.half_size,
            this.position_x + this.half_size,
            this.position_y + this.half_size,
            this.position_z + this.half_size
        );
    }

    calculate_signed_distance(x, y, z) {
        const pos_x = Math.abs(x - this.position_x) - this.half_size + this.corner_radius;
        const pos_y = Math.abs(y - this.position_y) - this.half_size + this.corner_radius;
        const pos_z = Math.abs(z - this.position_z) - this.half_size + this.corner_radius;
        const side = Math.min(Math.max(pos_x, pos_y, pos_z), 0.0);
        const dx = Math.max(pos_x, 0);
        const dy = Math.max(pos_y, 0);
        const dz = Math.max(pos_z, 0);
        return Math.sqrt(dx*dx + dy*dy + dz*dz) + side - this.corner_radius;
    }
}

export class SignedDistanceFieldCuboid extends SignedDistanceField {
    constructor(position = [0, 0, 0], size = [ 20, 20, 20 ], corner_radius = 0) {
        super();
        this.position_x    = position[0];
        this.position_y    = position[1];
        this.position_z    = position[2];
        this.half_size_x   = 0.5 * size[0];
        this.half_size_y   = 0.5 * size[1];
        this.half_size_z   = 0.5 * size[2];
        this.corner_radius = corner_radius;
        this.bounding_box  = new BoundingBox(
            this.position_x - this.half_size_x,
            this.position_y - this.half_size_y,
            this.position_z - this.half_size_z,
            this.position_x + this.half_size_x,
            this.position_y + this.half_size_y,
            this.position_z + this.half_size_z,
        );
    }

    calculate_signed_distance(x, y, z) {
        const pos_x = Math.abs(x - this.position_x) - this.half_size_x + this.corner_radius;
        const pos_y = Math.abs(y - this.position_y) - this.half_size_y + this.corner_radius;
        const pos_z = Math.abs(z - this.position_z) - this.half_size_z + this.corner_radius;
        const side = Math.min(Math.max(pos_x, pos_y, pos_z), 0.0);
        const dx = Math.max(pos_x, 0);
        const dy = Math.max(pos_y, 0);
        const dz = Math.max(pos_z, 0);
        return Math.sqrt(dx*dx + dy*dy + dz*dz) + side - this.corner_radius;
    }
}

export class SignedDistanceFieldUnion extends SignedDistanceField {
    constructor(... objects) {
        super();
        this.objects = objects;
        this.bounding_box = BoundingBox.enclose_distance_fields(this.objects);
    }

    calculate_signed_distance(x, y, z) {
        let distance = this.objects[0].calculate_signed_distance(x, y, z);

        for (let idx = 1; idx < this.objects.length; ++idx) {
            distance = Math.min(distance, this.objects[idx].calculate_signed_distance(x, y, z));
        }

        return distance;
    }
}

export class SignedDistanceFieldDifference extends SignedDistanceField {
    constructor(... objects) {
        super();
        this.objects = objects;
        this.bounding_box = this.objects[0].bounding_box.clone();
    }

    calculate_signed_distance(x, y, z) {
        let distance = this.objects[0].calculate_signed_distance(x, y, z);

        for (let idx = 1; idx < this.objects.length; ++idx) {
            distance = Math.max(distance, -this.objects[idx].calculate_signed_distance(x, y, z));
        }

        return distance;
    }
}

export class SignedDistanceFieldIntersection extends SignedDistanceField {
    constructor(... objects) {
        super();
        this.objects = objects;
        this.bounding_box = BoundingBox.intersecting_distance_fields(this.objects);
    }

    calculate_signed_distance(x, y, z) {
        let distance = this.objects[0].calculate_signed_distance(x, y, z);

        for (let idx = 1; idx < this.objects.length; ++idx) {
            distance = Math.max(distance, this.objects[idx].calculate_signed_distance(x, y, z));
        }

        return distance;
    }
}

export class SignedDistanceFieldUnionSmooth extends SignedDistanceField {
    constructor(smoothness, ... objects) {
        super();
        this.smoothness = smoothness;
        this.objects = objects;
        this.bounding_box = BoundingBox.enclose_distance_fields(this.objects);
    }
    
    calculate_signed_distance(x, y, z) {
        let distance0 = this.objects[0].calculate_signed_distance(x, y, z);

        for (let idx = 1; idx < this.objects.length; ++idx) {
            let distance1 = this.objects[idx].calculate_signed_distance(x, y, z);
            let h = clamp(0.5 + 0.5 * (distance1 - distance0) / this.smoothness, 0, 1);
            let m = distance1 + (distance0 - distance1) * h;
            distance0 = m - this.smoothness * h * (1 - h);
        }

        return distance0;
    }
}

export class SignedDistanceFieldDifferenceSmooth extends SignedDistanceField {
    constructor(smoothness, ... objects) {
        super();
        this.smoothness = smoothness;
        this.objects = objects;
        this.bounding_box = this.objects[0].bounding_box.clone();
    }
    
    calculate_signed_distance(x, y, z) {
        let distance0 = this.objects[0].calculate_signed_distance(x, y, z);

        for (let idx = 1; idx < this.objects.length; ++idx) {
            let distance1 = this.objects[idx].calculate_signed_distance(x, y, z);
            let h = clamp(0.5 - 0.5 * (distance1 + distance0) / this.smoothness, 0, 1);
            let m = distance0 + (-distance1 - distance0) * h;
            distance0 = m + this.smoothness * h * (1 - h);
        }

        return distance0;
    }
}

export class SignedDistanceFieldIntersectionSmooth extends SignedDistanceField {
    constructor(smoothness, ... objects) {
        super();
        this.smoothness = smoothness;
        this.objects = objects;
        this.bounding_box = BoundingBox.intersecting_distance_fields(this.objects);
    }
    
    calculate_signed_distance(x, y, z) {
        let distance0 = this.objects[0].calculate_signed_distance(x, y, z);

        for (let idx = 1; idx < this.objects.length; ++idx) {
            let distance1 = this.objects[idx].calculate_signed_distance(x, y, z);
            let h = clamp(0.5 - 0.5 * (distance1 - distance0) / this.smoothness, 0, 1);
            let m = distance1 + (distance0 - distance1) * h;
            distance0 = m + this.smoothness * h * (1 - h);
        }

        return distance0;
    }
}

export class SignedDistanceFieldTranslate extends SignedDistanceField {
    constructor(x, y, z, field) {
        super();
        this.translate_x = x;
        this.translate_y = y;
        this.translate_z = z;
        this.field       = field;
        this.bounding_box = new BoundingBox(
            field.bounding_box.min_x + x,
            field.bounding_box.min_y + y,
            field.bounding_box.min_z + z,
            field.bounding_box.max_x + x,
            field.bounding_box.max_y + y,
            field.bounding_box.max_z + z,
        );
    }

    calculate_signed_distance(x, y, z) {
        return this.field.calculate_signed_distance(x - this.translate_x,
                                                    y - this.translate_y,
                                                    z - this.translate_z);
    }
}

export class SignedDistanceFieldTransform extends SignedDistanceField {
    static temp_vector = new Vector3(0, 0, 0);
    constructor(field, matrix) {
        super();
        this.field = field;
        this.matrix = matrix.clone();
        this.inverse_matrix = matrix.inverse();

        this.bounding_box = this.calculate_bounding_box();
    }

    calculate_signed_distance(x, y, z) {
        SignedDistanceFieldTransform.temp_vector.set(x, y, z);
        SignedDistanceFieldTransform.temp_vector.transform(this.inverse_matrix);
        console.log("inverse_matrix:", this.inverse_matrix);
        console.log("temp_vector:", SignedDistanceFieldTransform.temp_vector);
        return this.field.calculate_signed_distance(
            SignedDistanceFieldTransform.temp_vector.x,
            SignedDistanceFieldTransform.temp_vector.y,
            SignedDistanceFieldTransform.temp_vector.z);
    }
}

export class SignedDistanceFieldOffset extends SignedDistanceField {
    constructor(object, offset) {
        super();
        this.offset = offset;
        this.object = object;
        this.bounding_box = this.object.bounding_box.clone().grow(this.offset);
    }

    calculate_signed_distance(x, y, z) {
        return this.object.calculate_signed_distance(x, y, z) - this.offset;
    }
}

export function cube({ position = [0, 0, 0], size = 10, corner_radius }) {
    return new SignedDistanceFieldCube(position, size, corner_radius);
}

export function cuboid({ position = [0, 0, 0], size = [ 20, 20, 20 ], corner_radius = 0 }) {
    return new SignedDistanceFieldCuboid(position, size, corner_radius);
}

export function sphere({ position = [0, 0, 0], radius = 5 }) {
    return new SignedDistanceFieldSphere(position, radius);
}

export function cylinder({ position = [ 0, 0, 0 ], height = 20, radius = 10 }) {
    return new SignedDistanceFieldCylinder(position, height, radius);
}
