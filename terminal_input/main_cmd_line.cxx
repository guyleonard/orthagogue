/*
 * Copyright 2012 Ole Kristian Ekseth (oekseth@gmail.com)
 *
 * This file is part of orthAgogue.
 *
 * orthAgogue is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * orthAgogue is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with orthAgogue. If not, see <http://www.gnu.org/licenses/>.
 */
#include "main_cmd_line.h"
#include "cmd_list.h"
#include "cmd_argument.h"
#include "../configure.h"
/**
   ---------------------------------------------------------------------------------------
   @Name: main.cxx -- This file is regarded as the launcher of the software, implying using the command parsing library devloped.
   @Author: Ole Kristian Ekseth (oekseth)
   @Last_Major_Update: 28.12.2011 by oekseth (cleanup)
   ---------------------------------------------------------------------------------------
**/
int main() {
#ifdef assert_code
  const bool print_info = true;
  cmd_list::assert_class(print_info);
  cmd_argument::assert_class(print_info);
#endif
}
