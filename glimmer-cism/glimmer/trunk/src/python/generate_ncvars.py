#! /usr/bin/env python

# Copyright 2004, Magnus Hagdorn
# 
# This file is part of glimmer.
# 
# PyGMT is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# 
# PyGMT is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with PyGMT; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

# python script used to generate source code files given a variable definition file

import ConfigParser, sys, time, string,re, os.path

NOATTRIB = ['name','dimensions','dimlen','data','factor','load','hot','type']
hotvars = []
dimensions = {}
module = {}

def dimid(name):
    return '%s_dimid'%name

def is_dimvar(var):
    """Return True if variable is associated with a dimension.

    this is assumed to be the case if no time dim is present
    """

    if len(string.split(var['dimensions'],',')) == 1 and 'dimlen' in var:
        return True
    else:
        return False
    
class Variables(dict):
    """Dictionary containing variable definitions."""

    def __init__(self,filename):
        """Initialise Variable class.

        filename: name of file containing variable definitions."""

        dict.__init__(self)

        # reading variable configuration file
        vars = ConfigParser.ConfigParser()
        vars.readfp(open(filename))

        for v in vars.sections():
            if v == 'VARSET':
                for (name, value) in vars.items(v):
                    module[name]=value
                continue
            vardef = {}
            vardef['name'] = v
            for (name, value) in vars.items(v):
                vardef[name] = value
            if 'hot' in vardef:
                if vardef['hot'].lower() in ['1','true','t']:
                    hotvars.append(v)
                    vardef['load'] = '1'
            if 'type' not in vardef:
                vardef['type'] = 'float'
            # handle dims
            for d in vardef['dimensions'].split(','):
                d=d.strip()
                if 'dimlen' in vardef:
                    dimensions[d] = vardef['dimlen']
                if d not in dimensions:
                    dimensions[d] = '-1'
            self.__setitem__(v,vardef)

    def keys(self):
        """Reorder standard keys alphabetically."""
        dk = []
        vk = []
        for v in dict.keys(self):
            if is_dimvar(self.__getitem__(v)):
                dk.append(v)
            else:
                vk.append(v)
        dk.sort()
        vk.sort()
        return dk+vk

class PrintVars:
    """Base class for printing variables."""
    canhandle = None
    comment = '!'

    def __init__(self,filename,outname=None):
        """Initialise.

        filename: name of file to be processed."""
        if os.path.basename(filename) != self.canhandle:
            raise NotImplementedError, 'Can only handle %s'%self.canhandle

        self.infile = open(filename,'r')
        if outname==None:
            self.stream = open(self.canhandle[:-3],'w')
        else:
            self.stream = open(outname,'w')

        self.handletoken = {'!GENVARS!' : self.print_var}

    def print_warning(self):
        """Write a warning message to stream"""

        self.stream.write("%s\n"%(80*self.comment))
        self.stream.write("%s WARNING: this file was automatically generated on\n%s %s\n%s from %s\n"%(self.comment,
                                                                                                       self.comment,time.strftime("%a, %d %b %Y %H:%M:%S +0000", time.gmtime()),
                                                                                                       self.comment, self.canhandle))
        self.stream.write("%s\n\n"%(80*self.comment))
        
    def print_var(self, var):
        """Template for writing single variable"""

        raise NotImplementedError, 'You should use one of the derived classes'

    def write(self,vars):
        """Merge file with definitions"""

        self.print_warning()
        for l in self.infile.readlines():
            for token in self.handletoken:
                if string.find(l,token) is not -1:
                    break
            if string.find(l,token) is not -1:
                for v in vars.keys():
                    self.handletoken[token](vars[v])
            else:
                self.stream.write("%s"%l)
        self.infile.close()
        self.stream.close()

class PrintDoc(PrintVars):
    """Process varlist.tex"""
    canhandle = 'varlist.tex.in'
    comment = '%'

    def __init__(self,filename):
        """Initialise.

        filename: name of file to be processed."""

        PrintVars.__init__(self,filename,'%s_varlist.tex'%module['name'])
        
    def print_var(self, var):
        """Write single variable block to stream for ncdf_params."""

        # skip variables associated with dimension 
        load = ''
        if 'load' in var:
            if var['load'].lower() in ['1','true','t']:
                load = '$^\\ast$'

        self.stream.write("\\texttt{%s}%s & %s & %s\\\\\n"%(var['name'],load,var['long_name'],
                                                            var['units'].replace('_','\_')))
        if 'standard_name' in var:
            self.stream.write("&CF name: \\texttt{%s}&\\\\\n"%(var['standard_name'].replace('_','\_')))
        self.stream.write("\\hline\n")

class PrintNC_template(PrintVars):
    """Process ncdf_template.f90.in"""
    canhandle = 'ncdf_template.f90.in'
    
    def __init__(self,filename):
        """Initialise.

        filename: name of file to be processed."""

        PrintVars.__init__(self,filename,'%s_io.f90'%module['name'])
        self.numvars = 0
        self.handletoken['!GENVAR_VARDEF!'] = self.print_vardef
        self.handletoken['!GENVAR_WRITE!'] = self.print_var_write
        self.handletoken['!GENVAR_READ!'] = self.print_var_read

    def write(self,vars):
        """Merge ncdf.F90.in with definitions."""

        numvars = 0
        for v in vars:
            if vars[v]['dimensions'] != v:
                numvars = numvars + 1

        self.thisvar = 1

        self.print_warning()
        for l in self.infile.readlines():
            for k in module.keys():
                l = l.replace(k.upper(),module[k])
            for token in self.handletoken:
                if string.find(l,token) is not -1:
                    break
            if string.find(l,token) is not -1:
                for v in vars.keys():
                    self.handletoken[token](vars[v])
            elif '!GENVAR_HOT!' in l:
                self.print_varhot()
            elif '!GENVAR_DIMS!' in l:
                self.print_dimensions()
            elif '!GENVAR_CHECKDIM!' in l:
                self.print_checkdims()
            else:
                self.stream.write("%s"%l)
        self.infile.close()
        self.stream.close()

    def print_varhot(self):
        """Create list of hotstart variables."""

        hotvar = ''
        for v in hotvars:
            hotvar = hotvar + ' %s '%v
        self.stream.write("  character(len=*),private,parameter :: hotvars = '%s'\n"%hotvar)
        
    def print_vardef(self,var):
        """Write single variable block to stream for ncdf_file."""

        dims = string.split(var['dimensions'],',')
        dims.reverse()
        dimstring = dimid(dims[0].strip())
        for d in dims[1:]:
            dimstring = '%s, %s'%(dimstring,dimid(d.strip()))
        
        self.stream.write("    !     %s -- %s\n"%(var['name'],var['long_name'])) # writing comment
        spaces = 0
        idstring = 'varid'
        if not is_dimvar(var):
            spaces=3
            self.stream.write("    pos = index(NCO%%vars,' %s ')\n"%var['name'])
            self.stream.write("    status = nf90_inq_varid(NCO%%id,'%s',varid)\n"%var['name'])
            self.stream.write("    if (pos.ne.0) then\n")
            self.stream.write("      NCO%%vars(pos+1:pos+%d) = '%s'\n"%(len(var['name']),len(var['name'])*' '))
            self.stream.write("    end if\n")
            self.stream.write("    if (pos.ne.0 .and. status.eq.nf90_enotvar) then\n")
        self.stream.write("%s    call write_log('Creating variable %s')\n"%(spaces*' ',var['name']))
        self.stream.write("%s    status = nf90_def_var(NCO%%id,'%s',NF90_%s,(/%s/),%s)\n"%(spaces*' ',
                                                                                           var['name'],
                                                                                           var['type'].upper(),
                                                                                           dimstring,
                                                                                           idstring
                                                                                           ))
        self.stream.write("%s    call nc_errorhandle(__FILE__,__LINE__,status)\n"%(spaces*' '))
        for attrib in var:
            if attrib not in NOATTRIB:
                self.stream.write("%s    status = nf90_put_att(NCO%%id, %s, '%s', '%s')\n"%(spaces*' ',
                                                                                            idstring,
                                                                                            attrib,
                                                                                            var[attrib]))
        if not is_dimvar(var):
            self.stream.write("%s    if (CFproj_allocated(model%%projection)) then\n"%(spaces*' '))
            self.stream.write("%s       status = nf90_put_att(NCO%%id, %s, 'grid_mapping',glimmer_nc_mapvarname)\n"%(spaces*' ',idstring))
            self.stream.write("%s    end if\n"%(spaces*' '))
            self.stream.write("%s  end if\n"%(spaces*' '))
        self.stream.write("\n")

    def print_dimensions(self):
        """Set up dimensions."""

        dims = dimensions.keys()
        dims.sort()
        # generate list of dimension ids
        for d in dims:
            self.stream.write("    integer :: %s_dimid\n"%d)
        # get dimension ids
        self.stream.write("\n    ! defining dimensions\n")
        for d in dims:
            if dimensions[d]!='-1': # create a new dimension
                self.stream.write("    status = nf90_def_dim(NCO%%id,'%s',%s,%s)\n"%(d,dimensions[d],dimid(d)))
            else:
                self.stream.write("    status = nf90_inq_dimid(NCO%%id,'%s',%s)\n"%(d,dimid(d)))
            self.stream.write("    call nc_errorhandle(__FILE__,__LINE__,status)\n")

    def print_checkdims(self):
        """Produce code for checking dimension sizes"""

        dims = dimensions.keys()
        dims.sort()
        for d in dims:
            if dimensions[d]!='-1':
                self.stream.write("    status = nf90_inq_dimid(NCI%%id,'%s',dimid)\n"%(d))
                self.stream.write("    if (dimid.gt.0) then\n")
                self.stream.write("       status = nf90_inquire_dimension(NCI%id,dimid,len=dimsize)\n")
                self.stream.write("       if (dimsize.ne.%s) then\n"%dimensions[d])
                self.stream.write("          write(message,*) 'Error, reading file ',trim(NCI%%filename),' size %s does not match: ',%s\n"%(d,dimensions[d]))
                self.stream.write("          call write_log(message,GM_FATAL)\n")
                self.stream.write("       end if\n")
                self.stream.write("    end if\n")
                
        
    def print_var_write(self,var):
        """Write single variable block to stream for ncdf_file."""

        # skip variables associated with dimension 
        if not is_dimvar(var):
            dims = string.split(var['dimensions'],',')
            dims.reverse()
            for i in range(0,len(dims)):
                dims[i] = dims[i].strip()
            self.stream.write("    status = nf90_inq_varid(NCO%%id,'%s',varid)\n"%var['name'])
            self.stream.write("    if (status .eq. nf90_noerr) then\n")
            
            dimstring = ''
            spaces = ''
            for i in range(0,len(dims)):
                if i > 0:
                    dimstring = dimstring + ','
                if dims[i] == 'time':
                    dimstring = dimstring + 'outfile%timecounter'
                elif dims[i] == 'level':
                    dimstring = dimstring + 'up'
                else:
                    dimstring = dimstring + '1'
                
            if  'level' in dims:
                # handle 3D fields
                spaces = ' '*3
                self.stream.write("       do up=1,NCO%nlevel\n")

                        
            if 'factor' in var:
                data = '(%s)*(%s)'%(var['factor'], var['data'])
            else:
                data = var['data']
            self.stream.write("%s       status = nf90_put_var(NCO%%id, varid, &\n%s            %s, (/%s/))\n"%(spaces,
                                                                                                               spaces,data, dimstring))
            self.stream.write("%s       call nc_errorhandle(__FILE__,__LINE__,status)\n"%(spaces))

            if  'level' in dims:
                self.stream.write("       end do\n")
            # remove self since it's not time dependent
            if 'time' not in dims:
                self.stream.write("       NCO%%do_var(%s) = .False.\n"%(var_type(var)))
                
            self.stream.write("    end if\n\n")

    def print_var_read(self,var):
        """Write single variable block to stream for reading netCDF data."""

        if 'load' in var and not is_dimvar(var):
            if var['load'].lower() in ['1','true','t']:
                dims = string.split(var['dimensions'],',')
                dims.reverse()
                for i in range(0,len(dims)):
                    dims[i] = dims[i].strip()
                self.stream.write("    status = nf90_inq_varid(NCI%%id,'%s',varid)\n"%var['name'])
                self.stream.write("    if (status .eq. nf90_noerr) then\n")
                self.stream.write("       call write_log('  Loading %s')\n"%var['name'])
                dimstring = ''
                spaces = ''
                for i in range(0,len(dims)):
                    if i > 0:
                        dimstring = dimstring + ','
                    if dims[i] == 'time':
                        dimstring = dimstring + 'infile%current_time'
                    elif dims[i] == 'level':
                        dimstring = dimstring + 'up'
                    else:
                        dimstring = dimstring + '1'

                if  'level' in dims:
                    # handle 3D fields
                    spaces = ' '*3
                    self.stream.write("       do up=1,NCI%nlevel\n")

                self.stream.write("%s       status = nf90_get_var(NCI%%id, varid, &\n%s            %s, (/%s/))\n"%(spaces,
                                                                                                               spaces,var['data'], dimstring))
                self.stream.write("%s       call nc_errorhandle(__FILE__,__LINE__,status)\n"%(spaces))
                if 'factor' in var:
                    self.stream.write("%s       if (scale) %s = %s/(%s)\n"%(spaces,var['data'],var['data'],var['factor']))

                if  'level' in dims:
                    self.stream.write("       end do\n")
                
            self.stream.write("    end if\n\n")

def usage():
    """Short help message."""

    print 'Usage generate_ncvars.py vardef [outfile0.in [,... [,outfileN.in]]]'
    print 'generate source code files given a variable definition file'
    print ''
    print 'vardef: file containing variable definition'
    print 'outfile.in: output template to be processed'
    print 'print variables if no templates are given'

HandleFile={}
HandleFile['ncdf_template.f90.in'] = PrintNC_template
HandleFile['varlist.tex.in'] = PrintDoc

if __name__ == '__main__':

    if len(sys.argv) < 2:
        usage()
        sys.exit(1)

    vars = Variables(sys.argv[1])

    if len(sys.argv) == 2:
        for v in vars.keys():
            print v
            for o in vars[v]:
                print '%s: %s'%(o, vars[v][o])
            print ''
            print module
    else:
        for f in sys.argv[2:]:
            base_f = os.path.basename(f)
            if base_f in HandleFile:
                handle = HandleFile[base_f](f)
                handle.write(vars)
            else:
                handle = PrintNCDF_IO(f)
                handle.write(vars)
                     
