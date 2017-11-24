################################################################################
# Introduction
#
# Determine the original name of any of the arguments defined in a function call
# to do any opreations on them.
#
# match the parent call of wherever this function is used. The output of
# match.call() wil be presented as a list of # which we can access all elements
# we are interested in. By specifying sys.function() and sys.call() with values
# -1 we match the function call of the first parent frame (in this case
# function) of this get.args().
# Credit goes to: http://stackoverflow.com/questions/17256834/getting-the-arguments-of-a-parent-function-in-r-with-names
#
################################################################################
# License
#
# Copyright (C) 2016 Steffen Ehrmann
#
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# This program is distributed in the hope that it will be useful , but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# this program. If not, see <http://www.gnu.org/licenses/>.
#
################################################################################
# function body
get.args <- function () {
  as.list(
    match.call(
      definition = sys.function( -1 ),
      call = sys.call( -1 )
      )
    )[-1]
}
