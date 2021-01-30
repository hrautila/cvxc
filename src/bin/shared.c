/*
 * Copyright by libcvxc authors. See AUTHORS file in this archive.
 *
 * This file is part of libcvxc library. It is free software,
 * distributed under the terms of GNU Lesser General Public License Version 3, or
 * any later version. See the COPYING file included in this archive.
 */

#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include <dlfcn.h>
#include <regex.h>
#include "solver.h"

int solver_is_shared_name(const char *name)
{
    regex_t preq;
    char *re = "\\." SO_EXT "$";

    if (regcomp(&preq, re, REG_NOSUB) != 0)
        return 0;

    int match = regexec(&preq, name, 0, (regmatch_t *)0, 0) == 0;
    regfree(&preq);
    return match;
}

char *solver_type_from_name(char *name)
{
    regex_t preq;
    regmatch_t group[2];
    const char *re = ".*(conelp|cpl|cp)$";

    if (regcomp(&preq, re, REG_EXTENDED) != 0) {
        return 0;
    }

    char *p = 0;
    if (regexec(&preq, name, 2, group, 0) == 0) {
        if (group[1].rm_so != -1) {
            p = &name[group[1].rm_so];
        }
    }
    regfree(&preq);
    return p;
}

char *solver_find_shared(const char *name, const char *var)
{
    size_t pathlen, namelen;
    char *paths = 0;
    if (var)
        paths = getenv(var);
    if (!paths)
        paths = getenv("CVXC_LIBRARY_PATH");
    if (!paths)
        paths = getenv("LD_LIBRARY_PATH");
    if (!paths) {
        return access(name, R_OK) == 0 ? strdup(name) : (char *)0;
    }

    pathlen = strlen(paths);
    namelen = strlen(name);

    char *buf = malloc(pathlen + namelen + 2);
    char *s = paths;
    char *p = 0;
    int has_extension = solver_is_shared_name(name);
    do {
        p = strchr(s, ':');
        if (!p)
            p = strchr(s, '\0');

        if (!p)
            break;

        pathlen = p - s;
        if (pathlen == 0) {
            // adjacent colons means current directory
            buf[0] = '.';
            buf[1] = '/';
            pathlen = 2;
        } else {
            memcpy(buf, s, pathlen);
        }
        if (buf[pathlen - 1] != '/') {
            // last character of path is not '/'
            buf[pathlen] = '/';
            pathlen++;
        }
        if (has_extension) {
            memcpy(&buf[pathlen], name, namelen);
            buf[pathlen + namelen] = '\0';
        } else {
            sprintf(&buf[pathlen], "%s.%s", name, SO_EXT);
        }
        if (access(buf, R_OK) == 0)
            return buf;

        s = p + 1;
    } while (*p);

    free(buf);
    return 0;
}

void *solver_load_shared(const char *name)
{
    if (!name)
        return 0;

    char *filepath = solver_find_shared(name, 0);
    if (!filepath)
        return 0;

    void *dlh = dlopen(filepath, RTLD_NOW);
    if (!dlh) {
        fprintf(stderr, "error: %s - %s\n", name, dlerror());
    }
    free(filepath);
    return dlh;
}
