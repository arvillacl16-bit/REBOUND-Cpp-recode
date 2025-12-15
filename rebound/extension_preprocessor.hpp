/*
 * This file is part of the C++ translation of REBOUND.
 *
 * Original REBOUND (C) code by Hanno Rein and others.
 * This translation is licensed under the GNU General Public License v3 or later.
 *
 * REBOUND is free software: you can redistribute it and/or modify it under the terms
 * of the GNU General Public License as published by the Free Software Foundation,
 * either version 3 of the License, or (at your option) any later version.
 *
 * REBOUND is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with this program.
 * If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef FRIEND_PP_HPP
#define FRIEND_PP_HPP

/* --------------------------------
   Minimal config
   -------------------------------- */
#define FRIEND_PP_VARIADICS 1

/* --------------------------------
   Helper macros
   -------------------------------- */

// Get first argument
#define FRIEND_PP_HEAD(x, ...) x
// Get rest of arguments
#define FRIEND_PP_TAIL(x, ...) __VA_ARGS__
// Check if there are more arguments
#define FRIEND_PP_HAS_ARGS(...) FRIEND_PP_HAS_ARGS_IMPL(__VA_ARGS__, 1,0)
#define FRIEND_PP_HAS_ARGS_IMPL(_x, _y, N, ...) N

/* --------------------------------
   Recursive mapping
   -------------------------------- */

// Expand macro for first argument
#define FRIEND_PP_MAP_1(M, x) M(x)

// Expand macro for variadic arguments recursively
#define FRIEND_PP_MAP(M, ...) \
    M(FRIEND_PP_HEAD(__VA_ARGS__)) \
    IF(FRIEND_PP_HAS_ARGS(FRIEND_PP_TAIL(__VA_ARGS__)))(FRIEND_PP_MAP(M, FRIEND_PP_TAIL(__VA_ARGS__)), )

/* Helper IF macro */
#define IF(cond) _IF_##cond
#define _IF_1(t, f) t
#define _IF_0(t, f) f

/* --------------------------------
   Friend macros
   -------------------------------- */

#define FRIEND_CLASS(x) friend class x;
#define FRIEND_FUNCTION(x) friend x;

/* --------------------------------
   User macros
   -------------------------------- */

#define ADD_FRIENDS_CLASS(...) FRIEND_PP_MAP(FRIEND_CLASS, __VA_ARGS__)
#define ADD_FRIENDS_FUNC(...)  FRIEND_PP_MAP(FRIEND_FUNCTION, __VA_ARGS__)

#endif // FRIEND_PP_HPP
