/*
 * minjson - Minimal JSON parser/builder for JSON-RPC 2.0
 * 
 * This is a lightweight JSON library specifically designed for the DBGPS3
 * server's JSON-RPC protocol. It supports:
 *   - Parsing JSON objects, arrays, strings, numbers, booleans, null
 *   - Building JSON values programmatically
 *   - Serializing to string
 *   - Getting/setting items by key in objects
 *   - Getting items by index in arrays
 *
 * License: MIT
 */

#ifndef MINJSON_H
#define MINJSON_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* JSON value types */
typedef enum {
    MJ_NULL = 0,
    MJ_BOOL,
    MJ_NUMBER,
    MJ_STRING,
    MJ_ARRAY,
    MJ_OBJECT
} mj_type;

/* Forward declaration */
typedef struct mj_val_s mj_val;

/* A JSON value */
struct mj_val_s {
    mj_type type;
    union {
        int boolean;
        double number;
        char *string;
        struct { mj_val *items; int count; int capacity; } array;
        struct { char **keys; mj_val *vals; int count; int capacity; } object;
    };
};

/* ---- Parsing ---- */

/**
 * Parse a JSON string into a value.
 * Returns a new value on success, or {type=MJ_NULL} on failure.
 * Caller must free with mj_free().
 */
mj_val mj_parse(const char *json);

/* ---- Building ---- */

/** Create a null value */
mj_val mj_null(void);
/** Create a boolean value */
mj_val mj_bool(int val);
/** Create a number value */
mj_val mj_number(double val);
/** Create a string value (copies the string) */
mj_val mj_string(const char *str);
/** Create an empty array */
mj_val mj_array(void);
/** Create an empty object */
mj_val mj_object(void);

/* ---- Modifying ---- */

/** Append an item to an array. Returns 0 on success. */
int mj_array_append(mj_val *arr, mj_val item);

/**
 * Set a key in an object. If key exists, replaces the value.
 * Both key and value are copied. Returns 0 on success.
 */
int mj_object_set(mj_val *obj, const char *key, mj_val item);

/* ---- Accessing ---- */

/** Get an item from an object by key. Returns {type=MJ_NULL} if not found. */
mj_val mj_object_get(const mj_val *obj, const char *key);

/** Check if an object has a key */
int mj_object_has(const mj_val *obj, const char *key);

/** Get an item from an array by index. Returns {type=MJ_NULL} if out of bounds. */
mj_val mj_array_get(const mj_val *arr, int index);

/* ---- Serializing ---- */

/**
 * Serialize a value to a JSON string.
 * Returns a malloc'd string that the caller must free.
 * Returns NULL on failure.
 */
char *mj_serialize(const mj_val *val);

/**
 * Serialize a value to a JSON string with pretty printing.
 * Returns a malloc'd string that the caller must free.
 */
char *mj_serialize_pretty(const mj_val *val);

/* ---- Cleanup ---- */

/**
 * Deep-free a value and all its children.
 * Safe to call on stack-allocated mj_val structs.
 */
void mj_free(mj_val *val);

/**
 * Get a string value, returning default if not a string.
 */
const char *mj_get_string(const mj_val *val, const char *def);

/**
 * Get a number value, returning default if not a number.
 */
double mj_get_number(const mj_val *val, double def);

/**
 * Get an integer value, returning default if not a number.
 */
int mj_get_int(const mj_val *val, int def);

#ifdef __cplusplus
}
#endif

#endif /* MINJSON_H */
