/*
 * minjson - Minimal JSON parser/builder implementation
 * License: MIT
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "minjson.h"

/* ---- Helpers ---- */

static char *mj_strdup(const char *s)
{
    if (!s) return NULL;
    size_t len = strlen(s);
    char *d = (char*)malloc(len + 1);
    if (d) memcpy(d, s, len + 1);
    return d;
}

/* ---- Constructors ---- */

mj_val mj_null(void) { mj_val v; memset(&v, 0, sizeof(v)); v.type = MJ_NULL; return v; }
mj_val mj_bool(int val) { mj_val v; memset(&v, 0, sizeof(v)); v.type = MJ_BOOL; v.boolean = val; return v; }
mj_val mj_number(double val) { mj_val v; memset(&v, 0, sizeof(v)); v.type = MJ_NUMBER; v.number = val; return v; }
mj_val mj_string(const char *str) { mj_val v; memset(&v, 0, sizeof(v)); v.type = MJ_STRING; v.string = mj_strdup(str); return v; }

mj_val mj_array(void)
{
    mj_val v;
    memset(&v, 0, sizeof(v));
    v.type = MJ_ARRAY;
    v.array.capacity = 8;
    v.array.items = (mj_val*)calloc(v.array.capacity, sizeof(mj_val));
    v.array.count = 0;
    return v;
}

mj_val mj_object(void)
{
    mj_val v;
    memset(&v, 0, sizeof(v));
    v.type = MJ_OBJECT;
    v.object.capacity = 8;
    v.object.keys = (char**)calloc(v.object.capacity, sizeof(char*));
    v.object.vals = (mj_val*)calloc(v.object.capacity, sizeof(mj_val));
    v.object.count = 0;
    return v;
}

/* ---- Modifying ---- */

int mj_array_append(mj_val *arr, mj_val item)
{
    if (!arr || arr->type != MJ_ARRAY) return -1;
    if (arr->array.count >= arr->array.capacity) {
        arr->array.capacity = arr->array.capacity * 2;
        arr->array.items = (mj_val*)realloc(arr->array.items, arr->array.capacity * sizeof(mj_val));
        if (!arr->array.items) return -1;
    }
    arr->array.items[arr->array.count++] = item;
    return 0;
}

int mj_object_set(mj_val *obj, const char *key, mj_val item)
{
    if (!obj || obj->type != MJ_OBJECT || !key) return -1;
    /* Check if key already exists */
    int i;
    for (i = 0; i < obj->object.count; i++) {
        if (strcmp(obj->object.keys[i], key) == 0) {
            mj_free(&obj->object.vals[i]);
            obj->object.vals[i] = item;
            return 0;
        }
    }
    /* Add new key */
    if (obj->object.count >= obj->object.capacity) {
        obj->object.capacity = obj->object.capacity * 2;
        obj->object.keys = (char**)realloc(obj->object.keys, obj->object.capacity * sizeof(char*));
        obj->object.vals = (mj_val*)realloc(obj->object.vals, obj->object.capacity * sizeof(mj_val));
        if (!obj->object.keys || !obj->object.vals) return -1;
    }
    obj->object.keys[obj->object.count] = mj_strdup(key);
    obj->object.vals[obj->object.count] = item;
    obj->object.count++;
    return 0;
}

/* ---- Accessing ---- */

mj_val mj_object_get(const mj_val *obj, const char *key)
{
    if (!obj || obj->type != MJ_OBJECT || !key) return mj_null();
    int i;
    for (i = 0; i < obj->object.count; i++) {
        if (strcmp(obj->object.keys[i], key) == 0)
            return obj->object.vals[i];
    }
    return mj_null();
}

int mj_object_has(const mj_val *obj, const char *key)
{
    if (!obj || obj->type != MJ_OBJECT || !key) return 0;
    int i;
    for (i = 0; i < obj->object.count; i++) {
        if (strcmp(obj->object.keys[i], key) == 0) return 1;
    }
    return 0;
}

mj_val mj_array_get(const mj_val *arr, int index)
{
    if (!arr || arr->type != MJ_ARRAY || index < 0 || index >= arr->array.count)
        return mj_null();
    return arr->array.items[index];
}

const char *mj_get_string(const mj_val *val, const char *def)
{
    return (val && val->type == MJ_STRING && val->string) ? val->string : def;
}

double mj_get_number(const mj_val *val, double def)
{
    return (val && val->type == MJ_NUMBER) ? val->number : def;
}

int mj_get_int(const mj_val *val, int def)
{
    return (val && val->type == MJ_NUMBER) ? (int)val->number : def;
}

/* ---- Cleanup ---- */

void mj_free(mj_val *val)
{
    if (!val) return;
    int i;
    switch (val->type) {
    case MJ_STRING:
        free(val->string);
        val->string = NULL;
        break;
    case MJ_ARRAY:
        for (i = 0; i < val->array.count; i++)
            mj_free(&val->array.items[i]);
        free(val->array.items);
        val->array.items = NULL;
        val->array.count = 0;
        break;
    case MJ_OBJECT:
        for (i = 0; i < val->object.count; i++) {
            free(val->object.keys[i]);
            mj_free(&val->object.vals[i]);
        }
        free(val->object.keys);
        free(val->object.vals);
        val->object.keys = NULL;
        val->object.vals = NULL;
        val->object.count = 0;
        break;
    default:
        break;
    }
    val->type = MJ_NULL;
}

/* ---- Parser ---- */

typedef struct {
    const char *json;
    int pos;
    int len;
} mj_parser;

static void mj_skip_whitespace(mj_parser *p)
{
    while (p->pos < p->len) {
        char c = p->json[p->pos];
        if (c == ' ' || c == '\t' || c == '\n' || c == '\r')
            p->pos++;
        else
            break;
    }
}

static mj_val mj_parse_value(mj_parser *p);

static char *mj_parse_string_raw(mj_parser *p)
{
    if (p->pos >= p->len || p->json[p->pos] != '"') return NULL;
    p->pos++; /* skip opening quote */

    /* First pass: compute length */
    int start = p->pos;
    int out_len = 0;
    int i = start;
    while (i < p->len && p->json[i] != '"') {
        if (p->json[i] == '\\') {
            i++; /* skip escape char */
            if (i < p->len) {
                switch (p->json[i]) {
                case '"': case '\\': case '/': case 'b': case 'f':
                case 'n': case 'r': case 't':
                    out_len++; break;
                case 'u':
                    out_len += 4; /* simplified: keep as \uXXXX */
                    i += 4;
                    break;
                default:
                    out_len++; break;
                }
            }
            i++;
        } else {
            out_len++;
            i++;
        }
    }

    /* Allocate and copy */
    char *result = (char*)malloc(out_len + 1);
    int j = 0;
    i = start;
    while (i < p->len && p->json[i] != '"') {
        if (p->json[i] == '\\') {
            i++;
            if (i < p->len) {
                switch (p->json[i]) {
                case '"':  result[j++] = '"'; break;
                case '\\': result[j++] = '\\'; break;
                case '/':  result[j++] = '/'; break;
                case 'b':  result[j++] = '\b'; break;
                case 'f':  result[j++] = '\f'; break;
                case 'n':  result[j++] = '\n'; break;
                case 'r':  result[j++] = '\r'; break;
                case 't':  result[j++] = '\t'; break;
                case 'u':
                    /* Simplified: copy \uXXXX as-is */
                    result[j++] = '\\'; result[j++] = 'u';
                    i++;
                    { int k; for(k=0; k<4 && i<p->len; k++, i++) result[j++] = p->json[i]; }
                    continue;
                default: result[j++] = p->json[i]; break;
                }
            }
            i++;
        } else {
            result[j++] = p->json[i++];
        }
    }
    result[j] = '\0';

    if (i < p->len && p->json[i] == '"') i++; /* skip closing quote */
    p->pos = i;
    return result;
}

static mj_val mj_parse_string_val(mj_parser *p)
{
    mj_val v;
    memset(&v, 0, sizeof(v));
    v.type = MJ_STRING;
    v.string = mj_parse_string_raw(p);
    if (!v.string) { v.type = MJ_NULL; }
    return v;
}

static mj_val mj_parse_number_val(mj_parser *p)
{
    mj_val v = mj_null();
    if (p->pos >= p->len) return v;

    char *end;
    double num = strtod(p->json + p->pos, &end);
    if (end == p->json + p->pos) return v; /* no conversion */
    v.type = MJ_NUMBER;
    v.number = num;
    p->pos = (int)(end - p->json);
    return v;
}

static mj_val mj_parse_array(mj_parser *p)
{
    mj_val arr = mj_array();
    if (p->pos >= p->len || p->json[p->pos] != '[') return arr;
    p->pos++; /* skip [ */

    mj_skip_whitespace(p);
    if (p->pos < p->len && p->json[p->pos] == ']') {
        p->pos++;
        return arr; /* empty array */
    }

    while (p->pos < p->len) {
        mj_skip_whitespace(p);
        mj_val item = mj_parse_value(p);
        mj_array_append(&arr, item);
        mj_skip_whitespace(p);
        if (p->pos < p->len && p->json[p->pos] == ',') {
            p->pos++;
            continue;
        }
        break;
    }

    mj_skip_whitespace(p);
    if (p->pos < p->len && p->json[p->pos] == ']') p->pos++;
    return arr;
}

static mj_val mj_parse_object(mj_parser *p)
{
    mj_val obj = mj_object();
    if (p->pos >= p->len || p->json[p->pos] != '{') return obj;
    p->pos++; /* skip { */

    mj_skip_whitespace(p);
    if (p->pos < p->len && p->json[p->pos] == '}') {
        p->pos++;
        return obj; /* empty object */
    }

    while (p->pos < p->len) {
        mj_skip_whitespace(p);
        if (p->pos >= p->len || p->json[p->pos] != '"') break;

        char *key = mj_parse_string_raw(p);
        if (!key) break;

        mj_skip_whitespace(p);
        if (p->pos < p->len && p->json[p->pos] == ':') p->pos++; /* skip : */
        mj_skip_whitespace(p);

        mj_val val = mj_parse_value(p);
        mj_object_set(&obj, key, val);
        free(key);

        mj_skip_whitespace(p);
        if (p->pos < p->len && p->json[p->pos] == ',') {
            p->pos++;
            continue;
        }
        break;
    }

    mj_skip_whitespace(p);
    if (p->pos < p->len && p->json[p->pos] == '}') p->pos++;
    return obj;
}

static mj_val mj_parse_value(mj_parser *p)
{
    mj_skip_whitespace(p);
    if (p->pos >= p->len) return mj_null();

    char c = p->json[p->pos];

    if (c == '"') return mj_parse_string_val(p);
    if (c == '{') return mj_parse_object(p);
    if (c == '[') return mj_parse_array(p);
    if (c == 't') {
        if (p->pos + 3 < p->len && strncmp(p->json + p->pos, "true", 4) == 0) {
            p->pos += 4;
            return mj_bool(1);
        }
    }
    if (c == 'f') {
        if (p->pos + 4 < p->len && strncmp(p->json + p->pos, "false", 5) == 0) {
            p->pos += 5;
            return mj_bool(0);
        }
    }
    if (c == 'n') {
        if (p->pos + 3 < p->len && strncmp(p->json + p->pos, "null", 4) == 0) {
            p->pos += 4;
            return mj_null();
        }
    }
    if (c == '-' || (c >= '0' && c <= '9')) {
        return mj_parse_number_val(p);
    }

    return mj_null();
}

mj_val mj_parse(const char *json)
{
    if (!json) return mj_null();
    mj_parser p;
    p.json = json;
    p.pos = 0;
    p.len = (int)strlen(json);
    return mj_parse_value(&p);
}

/* ---- Serializer ---- */

typedef struct {
    char *buf;
    int len;
    int capacity;
} mj_writer;

static void mw_init(mj_writer *w)
{
    w->capacity = 256;
    w->buf = (char*)malloc(w->capacity);
    w->len = 0;
    w->buf[0] = '\0';
}

static void mw_ensure(mj_writer *w, int need)
{
    if (w->len + need >= w->capacity) {
        while (w->len + need >= w->capacity)
            w->capacity *= 2;
        w->buf = (char*)realloc(w->buf, w->capacity);
    }
}

static void mw_putc(mj_writer *w, char c)
{
    mw_ensure(w, 2);
    w->buf[w->len++] = c;
    w->buf[w->len] = '\0';
}

static void mw_puts(mj_writer *w, const char *s)
{
    int slen = (int)strlen(s);
    mw_ensure(w, slen + 1);
    memcpy(w->buf + w->len, s, slen);
    w->len += slen;
    w->buf[w->len] = '\0';
}

static void mw_write_string(mj_writer *w, const char *s)
{
    mw_putc(w, '"');
    if (s) {
        int i;
        for (i = 0; s[i]; i++) {
            switch (s[i]) {
            case '"':  mw_puts(w, "\\\""); break;
            case '\\': mw_puts(w, "\\\\"); break;
            case '\b': mw_puts(w, "\\b"); break;
            case '\f': mw_puts(w, "\\f"); break;
            case '\n': mw_puts(w, "\\n"); break;
            case '\r': mw_puts(w, "\\r"); break;
            case '\t': mw_puts(w, "\\t"); break;
            default:
                if ((unsigned char)s[i] < 0x20) {
                    char esc[8];
                    snprintf(esc, sizeof(esc), "\\u%04x", (unsigned char)s[i]);
                    mw_puts(w, esc);
                } else {
                    mw_putc(w, s[i]);
                }
                break;
            }
        }
    }
    mw_putc(w, '"');
}

static void mj_serialize_val(const mj_val *val, mj_writer *w);

static void mj_serialize_object(const mj_val *val, mj_writer *w, int indent, int pretty)
{
    mw_putc(w, '{');
    int i;
    for (i = 0; i < val->object.count; i++) {
        if (i > 0) { mw_putc(w, ','); }
        if (pretty) {
            mw_putc(w, '\n');
            int j; for(j=0; j<=indent; j++) mw_putc(w, ' ');
        }
        mw_write_string(w, val->object.keys[i]);
        mw_putc(w, ':');
        if (pretty) mw_putc(w, ' ');
        mj_serialize_val(&val->object.vals[i], w);
    }
    if (pretty && val->object.count > 0) {
        mw_putc(w, '\n');
        { int j; for(j=0; j<indent; j++) mw_putc(w, ' '); }
    }
    mw_putc(w, '}');
}

static void mj_serialize_array(const mj_val *val, mj_writer *w, int indent, int pretty)
{
    mw_putc(w, '[');
    int i;
    for (i = 0; i < val->array.count; i++) {
        if (i > 0) { mw_putc(w, ','); }
        if (pretty) { mw_putc(w, '\n'); int j; for(j=0; j<=indent; j++) mw_putc(w, ' '); }
        mj_serialize_val(&val->array.items[i], w);
    }
    if (pretty && val->array.count > 0) {
        mw_putc(w, '\n');
        { int j; for(j=0; j<indent; j++) mw_putc(w, ' '); }
    }
    mw_putc(w, ']');
}

static void mj_serialize_val(const mj_val *val, mj_writer *w)
{
    if (!val) { mw_puts(w, "null"); return; }
    switch (val->type) {
    case MJ_NULL:
        mw_puts(w, "null");
        break;
    case MJ_BOOL:
        mw_puts(w, val->boolean ? "true" : "false");
        break;
    case MJ_NUMBER: {
        char numbuf[64];
        /* Print integer-valued doubles without decimal */
        if (val->number == (double)(long long)val->number && !isinf(val->number) && !isnan(val->number)) {
            snprintf(numbuf, sizeof(numbuf), "%lld", (long long)val->number);
        } else {
            snprintf(numbuf, sizeof(numbuf), "%.6g", val->number);
        }
        mw_puts(w, numbuf);
        break;
    }
    case MJ_STRING:
        mw_write_string(w, val->string);
        break;
    case MJ_ARRAY:
        mj_serialize_array(val, w, 0, 0);
        break;
    case MJ_OBJECT:
        mj_serialize_object(val, w, 0, 0);
        break;
    }
}

static void mj_serialize_val_pretty(const mj_val *val, mj_writer *w, int indent)
{
    if (!val) { mw_puts(w, "null"); return; }
    switch (val->type) {
    case MJ_NULL:
        mw_puts(w, "null");
        break;
    case MJ_BOOL:
        mw_puts(w, val->boolean ? "true" : "false");
        break;
    case MJ_NUMBER: {
        char numbuf[64];
        if (val->number == (double)(long long)val->number && !isinf(val->number) && !isnan(val->number)) {
            snprintf(numbuf, sizeof(numbuf), "%lld", (long long)val->number);
        } else {
            snprintf(numbuf, sizeof(numbuf), "%.6g", val->number);
        }
        mw_puts(w, numbuf);
        break;
    }
    case MJ_STRING:
        mw_write_string(w, val->string);
        break;
    case MJ_ARRAY:
        mj_serialize_array(val, w, indent + 2, 1);
        break;
    case MJ_OBJECT:
        mj_serialize_object(val, w, indent + 2, 1);
        break;
    }
}

char *mj_serialize(const mj_val *val)
{
    mj_writer w;
    mw_init(&w);
    mj_serialize_val(val, &w);
    return w.buf;
}

char *mj_serialize_pretty(const mj_val *val)
{
    mj_writer w;
    mw_init(&w);
    mj_serialize_val_pretty(val, &w, 0);
    return w.buf;
}
