/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.uam.es)
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or   
 * (at your option) any later version.                                 
 *                                                                     
 * This program is distributed in the hope that it will be useful,     
 * but WITHOUT ANY WARRANTY; without even the implied warranty of      
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       
 * GNU General Public License for more details.                        
 *                                                                     
 * You should have received a copy of the GNU General Public License   
 * along with this program; if not, write to the Free Software         
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA            
 * 02111-1307  USA                                                     
 *                                                                     
 *  All comments concerning this program package may be sent to the    
 *  e-mail address 'xmipp@cnb.uam.es'                                  
 ***************************************************************************/

#include <XmippData/xmippProgs.hh>
#include <XmippData/xmippArgs.hh>
#include <Reconstruction/Programs/Prog_CorrectPhase.hh>

class Prog_CorrectPhase_Params: public Prog_parameters {
public:
   CorrectPhase_Params cpprm;
   FileName            out_ctf;
public:
   void read(int argc, char **argv) _THROW {
      Prog_parameters::read(argc,argv);
      cpprm.read(argc,argv);
      cpprm.produce_side_info();
      out_ctf=get_param(argc,argv,"-output_ctf","");
      allow_time_bar=FALSE;
   }
   
   void show() {
      Prog_parameters::show();
      cpprm.show();
   }

   void usage() {
      Prog_parameters::usage();
      cpprm.usage();
      cerr << "  [-output_ctf <Xmipp Fourier file>]: If not given, the \n"
           << "                              input one is rewritten\n"
	   << "                              If a selfile is used this\n"
	   << "                              parameter is used as an output root\n"
      ;
   }
};

bool process_img(ImageXmipp &img, const Prog_parameters *prm) {
   Prog_CorrectPhase_Params *eprm=(Prog_CorrectPhase_Params *) prm;
   matrix2D< complex<double> > fft;
   FileName       fn_ctf;
   if (eprm->cpprm.multiple_CTFs) {
      fn_ctf=eprm->cpprm.CTF_filename(img.name());
      if (img.name().get_number() != fn_ctf.get_number()) {
         cerr << "Cannot find CTF for image " << img.name() << endl;
         return FALSE;
      }
      eprm->cpprm.ctf.read_mask(fn_ctf);
      cerr << "Correcting " << img.name() << " with " << fn_ctf << endl;
   }
   FourierTransform(img(),fft);
   eprm->cpprm.correct(fft);
   InverseFourierTransform(fft,img());

   // Correct the CTF itself
   if (eprm->cpprm.multiple_CTFs) {
      FileName fn_ctf_out;
      if (eprm->out_ctf=="") fn_ctf_out=fn_ctf;
      else fn_ctf_out=eprm->out_ctf+ItoA(fn_ctf.get_number(),5)+".fft";

      eprm->cpprm.correct(eprm->cpprm.ctf.mask2D);
      eprm->cpprm.ctf.write_mask(fn_ctf_out,2);
   }
   return TRUE;
}

bool process_vol(VolumeXmipp &vol, const Prog_parameters *prm) {
   cerr << "This process is not intended for volumes\n";
   return FALSE;
}

int main (int argc, char **argv) {
   Prog_CorrectPhase_Params prm;
   SF_main(argc, argv, &prm, (void*)&process_img, (void*)&process_vol);
   if (!prm.cpprm.multiple_CTFs) {
      // Correct the CTF itself
      prm.cpprm.correct(prm.cpprm.ctf.mask2D);
      if (prm.out_ctf=="") prm.cpprm.ctf.write_mask(prm.cpprm.fn_ctf,2);
      else                 prm.cpprm.ctf.write_mask(prm.out_ctf,2);
   }
}
