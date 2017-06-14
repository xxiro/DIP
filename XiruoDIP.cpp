
//************************************************************************
//************************************************************************
// BertImage
// created: february 2006
// updated: may 2011
// How to load, display, process and save images with wxWindows
// Author: Pascal Bertolino, GIPSA-lab laboratory, Grenoble - France
// Email pascal.bertolino@gipsa-lab.fr
// Web http://www.gipsa-lab.inpg.fr/~pascal.bertolino/
// tested with xWidget 2.8.7 under Linux and Windows
//************************************************************************
//************************************************************************

// For compilers that support precompilation, includes "wx/wx.h".
#include "wx/wxprec.h"

#ifdef __BORLANDC__
    #pragma hdrstop
#endif

#ifndef WX_PRECOMP
    #include "wx/wx.h"
#endif

#include <wx/image.h>
#include <wx/file.h>
#include <wx/bitmap.h>
#include <wx/textdlg.h> 
#include <cstdlib>
#include <iostream>
#include <vector>
#include <math.h>
#include <cmath>
#define APP_NAME "Xiruo's Image Processing App"

using namespace std;

enum
{
    ID_QUIT = 1,
    ID_ABOUT,
	ID_LOAD,
	ID_SAVE,
	ID_MKGREY,
	ID_PROCESS,
	ID_INVERT,
	ID_RESCALE,
	ID_SHIFT,
	ID_AVERAGE,
	ID_LAPLACIAN,
	ID_NOISEINTRO,
	ID_MIN,
	ID_MAX,
	ID_MEDIAN,
	ID_NEGLIN,
	ID_LOG,
	ID_POW,
	ID_HIST,
	ID_AUTHRES,
	ID_RESET,
	ID_BEST_SIZE
};

//************************************************************************
//************************************************************************
// Canvas class (where we display the image)
//************************************************************************
//************************************************************************

//------------------------------------------------------------------------
class MyCanvas: public wxPanel
//------------------------------------------------------------------------
{
public:
    MyCanvas( wxWindow *parent, wxWindowID, const wxPoint &pos, const wxSize &size ) ;
    ~MyCanvas() ;
	void LoadImage(wxString fileName) ;
	void SaveImage(wxString fileName) ;
	void MkGrey();
	void InvertImage() ;
	void Rescale() ;
	void Shift(double value);
	void Average();
	void Laplacian();
	void Noise_intro(long num);
	void Min();
	void Max();
	void Median();
	void NegLin();
	void Log();
	void Pow();
	void Hist();
	void AuThres();

	void Reset();

	void BestSize() ;

private:
	int m_imageWidth ;
	int m_imageHeight ;
    wxBitmap m_imageBitmap ;	// used to display the image
	wxImage *m_imageRGB ;		// used to load the image
	wxImage *temp_imageRGB ;	// used as a buffer for convolutional filters

	unsigned char* m_imageData; 	// used to store the initial image data for reset
	unsigned char* temp_imageData;	// used as a buffer for convolutional filters
    void OnPaint(wxPaintEvent &event) ;

    DECLARE_EVENT_TABLE()
};

BEGIN_EVENT_TABLE(MyCanvas, wxPanel)
	EVT_PAINT(MyCanvas::OnPaint)
END_EVENT_TABLE()

//------------------------------------------------------------------------
MyCanvas::MyCanvas( wxWindow *parent, wxWindowID id,
                    const wxPoint &pos, const wxSize &size )
        : wxPanel( parent, id, pos, size, wxSUNKEN_BORDER)
//------------------------------------------------------------------------
{
	m_imageRGB = NULL ;
	temp_imageRGB = NULL;
}

//------------------------------------------------------------------------
MyCanvas::~MyCanvas()
//------------------------------------------------------------------------
{
	if (m_imageData)
		free(m_imageData) ;
	if (m_imageRGB)
		delete m_imageRGB ;
}

//------------------------------------------------------------------------
void MyCanvas::LoadImage(wxString fileName)
//------------------------------------------------------------------------
{
	if (m_imageData)
		free (m_imageData) ;
	if (m_imageRGB){
		delete m_imageRGB ;
		delete m_imageData ;
	}
// open image dialog box
	m_imageRGB = new wxImage(fileName, wxBITMAP_TYPE_ANY, -1); // ANY => can load many image formats
	m_imageBitmap = wxBitmap(*m_imageRGB, -1); // ...to get the corresponding bitmap

	m_imageWidth = m_imageRGB->GetWidth() ;
	m_imageHeight = m_imageRGB->GetHeight() ;

	m_imageData = ( unsigned char* ) malloc( m_imageWidth * m_imageHeight * 3 );
	memcpy(m_imageData, m_imageRGB->GetData(),m_imageWidth * m_imageHeight * 3);


// update GUI size
	SetSize(m_imageWidth, m_imageHeight) ;
	GetParent()->SetClientSize(GetSize()) ;

// update display
	Refresh(false) ;
}

//------------------------------------------------------------------------
void MyCanvas::SaveImage(wxString fileName)
//------------------------------------------------------------------------
{
	bool b ;

	//wxImage* tempImage = new wxImage(m_imageWidth, m_imageHeight, m_myImage, true); // lend my image buffer...
	//b = tempImage->SaveFile(fileName) ;
	//delete(tempImage) ;		// buffer not needed any more
	b = m_imageRGB->SaveFile(fileName) ;
	if(!b)
		wxMessageBox(wxT("A problem occured during saving"));
}


//------------------------------------------------------------------------


void MyCanvas::InvertImage() //simple inversion of image
{
    for( int i=0;i<m_imageWidth;i++)
       for(int j=0;j<m_imageHeight;j++){
 	m_imageRGB->SetRGB(i,j,255-m_imageRGB->GetRed(i,j), 
				255-m_imageRGB->GetGreen(i,j),
				255-m_imageRGB->GetBlue(i,j));
    }	
 	Refresh(false);
}

void MyCanvas::MkGrey() //make a color image into grey image
{
	double grey;

    for( int i=0;i<m_imageWidth;i++)
       for(int j=0;j<m_imageHeight;j++){
       		grey = 0.2989*m_imageRGB->GetRed(i,j)+
       			   0.5870*m_imageRGB->GetGreen(i,j)+
       			   0.1140*m_imageRGB->GetBlue(i,j);
 			m_imageRGB->SetRGB(i,j,grey,grey,grey);
    }	
 	Refresh(false);
}

void MyCanvas::Rescale() //scale up the brightness by a factor of 2.5
{

    for( int i=0;i<m_imageWidth;i++)
       for(int j=0;j<m_imageHeight;j++){

       		int r = (int)(2.5*(m_imageRGB->GetRed(i,j)+0.5));
       		int g = (int)(2.5*(m_imageRGB->GetGreen(i,j)+0.5));
       		int b = (int)(2.5*(m_imageRGB->GetBlue(i,j)+0.5));

 			m_imageRGB->SetRGB(i,j,r < 255 ? r : 255, 
								   g < 255 ? g : 255,
								   b < 255 ? b : 255);
    }

 	Refresh(false);    	
}

void MyCanvas::Shift(double value) //shift the brightness by an input amount
{

    for( int i=0;i<m_imageWidth;i++)
       for(int j=0;j<m_imageHeight;j++){

       		int r = (int)((m_imageRGB->GetRed(i,j)+value)+0.5);
       		int g = (int)((m_imageRGB->GetGreen(i,j)+value)+0.5);
       		int b = (int)((m_imageRGB->GetBlue(i,j)+value)+0.5);

 			m_imageRGB->SetRGB(i,j,(r > 255) ? 255 : ((r < 0) ? 0 : r), 
								   (g > 255) ? 255 : ((g < 0) ? 0 : g),
								   (b > 255) ? 255 : ((b < 0) ? 0 : b));
    }

 	Refresh(false);    	
}
//------------------------------------------------------------------------

void MyCanvas::Average() //average filter
{
	double average[3*3] = {
		1., 1., 1.,
		1., 1., 1.,
		1., 1., 1.,
	};

	int ix, iy, l;
	int kx, ky;
	double cp[3];
	double divisor = 9.;
	double offset = 0.;
	int Ks = 1;

	temp_imageData = (unsigned char*) malloc(m_imageWidth*m_imageHeight*3);
	memcpy(temp_imageData,m_imageRGB->GetData(),m_imageWidth*m_imageHeight*3);
	temp_imageRGB = new wxImage(m_imageWidth, m_imageHeight, temp_imageData, false);

    for(ix=Ks; ix < m_imageWidth -Ks; ix++) {
    	for(iy=Ks; iy< m_imageHeight-Ks; iy++) {
			cp[0] = cp[1] = cp[2] = 0.0;
			for(kx=-Ks; kx <= Ks; kx++) {
	 			for(ky=-Ks; ky <= Ks; ky++) {

	 				cp[0] += (average[(kx+Ks) + (ky+Ks)*(2*Ks+1)]/divisor) * ((temp_imageRGB->GetRed(ix+kx,iy+ky)+offset));
	 				cp[1] += (average[(kx+Ks) + (ky+Ks)*(2*Ks+1)]/divisor) * ((temp_imageRGB->GetGreen(ix+kx,iy+ky)+offset));
	 				cp[2] += (average[(kx+Ks) + (ky+Ks)*(2*Ks+1)]/divisor) * ((temp_imageRGB->GetBlue(ix+kx,iy+ky)+offset)); 	 				
	  			}
			}		
			for(l=0; l<3; l++)
				cp[l] = (cp[l]>255.0) ? 255.0 : ((cp[l]<0.0) ? 0.0 : cp[l]) ;
		
			m_imageRGB->SetRGB(ix,iy,cp[0],cp[1],cp[2]);
      }
    }
	delete temp_imageRGB;
    Refresh(); 

}

void MyCanvas::Laplacian() //8-neighbor Laplacian filter
{
	double mat[3*3] = {
		-1., -1., -1.,
		-1., 8., -1.,
		-1., -1., -1.,
	};

	int ix, iy, l;
	int kx, ky;
	double cp[3];
	double divisor = 1.;
	double offset = 0.;
	int Ks = 1;

	temp_imageData = (unsigned char*) malloc(m_imageWidth*m_imageHeight*3);
	memcpy(temp_imageData,m_imageRGB->GetData(),m_imageWidth*m_imageHeight*3);
	temp_imageRGB = new wxImage(m_imageWidth, m_imageHeight, temp_imageData, false);

    for(ix=Ks; ix < m_imageWidth -Ks; ix++) {
    	for(iy=Ks; iy< m_imageHeight-Ks; iy++) {
			cp[0] = cp[1] = cp[2] = 0.0;
			for(kx=-Ks; kx <= Ks; kx++) {
	 			for(ky=-Ks; ky <= Ks; ky++) {

	 				cp[0] += (mat[(kx+Ks) + (ky+Ks)*(2*Ks+1)]/divisor) * ((temp_imageRGB->GetRed(ix+kx,iy+ky)+offset));
	 				cp[1] += (mat[(kx+Ks) + (ky+Ks)*(2*Ks+1)]/divisor) * ((temp_imageRGB->GetGreen(ix+kx,iy+ky)+offset));
	 				cp[2] += (mat[(kx+Ks) + (ky+Ks)*(2*Ks+1)]/divisor) * ((temp_imageRGB->GetBlue(ix+kx,iy+ky)+offset)); 	 				

	  			}
			}		
			for(l=0; l<3; l++)
				cp[l] = (cp[l]>255.0) ? 255.0 : ((cp[l]<0.0) ? 0.0 : cp[l]) ;
		
			m_imageRGB->SetRGB(ix,iy,cp[0],cp[1],cp[2]);
      }
    }
	delete temp_imageRGB;
    Refresh(); 

}

void MyCanvas::Noise_intro(long num) //introducing salt and pepper noise with input number of random iterations
{
	int ii, rx, ry;

	for (ii = 0; ii < num; ii++)
	{
		rx = rand() % m_imageWidth;
		ry = rand() % m_imageHeight;
		m_imageRGB->SetRGB(rx, ry, 255, 255, 255);
	};

	for (ii = 0; ii < num; ii++)
	{
		rx = rand() % m_imageWidth;
		ry = rand() % m_imageHeight;
		m_imageRGB->SetRGB(rx, ry, 0, 0, 0);
	};

	Refresh();	
}

void MyCanvas::Min() //minimum filter
{

	int ix, iy, kx, ky;
	int sx = 0;
	int sy = 0;
	int temp, smallest;
	int Ks = 1;
	bool pflag;
	
	temp_imageData = (unsigned char*) malloc(m_imageWidth*m_imageHeight*3);
	memcpy(temp_imageData,m_imageRGB->GetData(),m_imageWidth*m_imageHeight*3);
	temp_imageRGB = new wxImage(m_imageWidth, m_imageHeight, temp_imageData, false);	

    for(ix=0; ix < m_imageWidth; ix++) {
    	for(iy=0; iy< m_imageHeight; iy++) {
    		smallest = temp_imageRGB->GetRed(ix,iy);
			for(kx=-Ks; kx <= Ks; kx++) {
	 			for(ky=-Ks; ky <= Ks; ky++) {

	 				pflag = ((ix+kx)>=0 && (iy+ky)>=0 && (ix+kx)<m_imageWidth && (iy+ky)<m_imageHeight); //

	 				if (pflag){
	 					temp = temp_imageRGB->GetRed(ix+kx,iy+ky);
	 					if(smallest>temp)
	 						smallest = temp;
	 				};
	  			};
			};
		
			m_imageRGB->SetRGB(ix,iy,smallest,smallest,smallest);
    	};
    };

	delete temp_imageRGB;

    Refresh();

}

void MyCanvas::Max() //maximum filter
{

	int ix, iy, kx, ky;
	int sx = 0;
	int sy = 0;
	int temp, biggest;
	int Ks = 1;
	
	temp_imageData = (unsigned char*) malloc(m_imageWidth*m_imageHeight*3);
	memcpy(temp_imageData,m_imageRGB->GetData(),m_imageWidth*m_imageHeight*3);
	temp_imageRGB = new wxImage(m_imageWidth, m_imageHeight, temp_imageData, false);	

    for(ix=Ks; ix < (m_imageWidth-Ks); ix++) {
    	for(iy=Ks; iy< (m_imageHeight-Ks); iy++) {
    		biggest = temp_imageRGB->GetRed(ix,iy);
			for(kx=-Ks; kx <= Ks; kx++) {
	 			for(ky=-Ks; ky <= Ks; ky++) {

	 					temp = temp_imageRGB->GetRed(ix+kx,iy+ky);
	 					if(biggest<temp)
	 						biggest = temp;
	  			};
			};
		
			m_imageRGB->SetRGB(ix,iy,biggest,biggest,biggest);
    	};
    };

	delete temp_imageRGB;

    Refresh();

}

void MyCanvas::Median() //median filter
{

	int ix, iy, kx, ky;
	int sx = 0;
	int sy = 0;
	int temp, biggest, ii;
	int Ks = 1;

	vector<double> vec;

	temp_imageData = (unsigned char*) malloc(m_imageWidth*m_imageHeight*3);
	memcpy(temp_imageData,m_imageRGB->GetData(),m_imageWidth*m_imageHeight*3);
	temp_imageRGB = new wxImage(m_imageWidth, m_imageHeight, temp_imageData, false);	

    for(ix=Ks; ix < (m_imageWidth-Ks); ix++) {
    	for(iy=Ks; iy< (m_imageHeight-Ks); iy++) {
    		
			for(kx=-Ks; kx <= Ks; kx++) {
	 			for(ky=-Ks; ky <= Ks; ky++) {

	 					vec.push_back(temp_imageRGB->GetRed(ix+kx,iy+ky));
	 					
	  			};
			};
			nth_element(vec.begin(), vec.begin()+vec.size()/2, vec.end());
			m_imageRGB->SetRGB(ix,iy,vec[vec.size()/2],vec[vec.size()/2],vec[vec.size()/2]);

			vec.clear();
    	};
    };

	delete temp_imageRGB;

    Refresh();

}

void MyCanvas::NegLin() //negative linear process (inversion here)
{

	int L = 256;
	int ix, iy;
	double s;
	for (ix = 0; ix< m_imageWidth; ix++)
		for (iy = 0; iy< m_imageHeight; iy++)
		{
			s = L - 1 - m_imageRGB->GetRed(ix,iy);
			m_imageRGB->SetRGB(ix,iy,s,s,s);
		}
			
	Refresh();


}

void MyCanvas::Log() //logarithm process
{

	int ix, iy;
	double s;
	double max = 0;
	for (ix = 0; ix< m_imageWidth; ix++)
		for (iy = 0; iy< m_imageHeight; iy++)
		{
			if (max<m_imageRGB->GetRed(ix,iy))
				max = m_imageRGB->GetRed(ix,iy);
		}
			

	for (ix = 0; ix< m_imageWidth; ix++)
		for (iy = 0; iy< m_imageHeight; iy++)
		{
			s = 255.0/log(max)*log(1.0+ m_imageRGB->GetRed(ix,iy));
			m_imageRGB->SetRGB(ix,iy,s,s,s);
		}
			
	Refresh();


}

void MyCanvas::Pow() //power law process (gamma)
{

	double p = 2;
	int ix, iy;
	double s;
	double max = 0;
	for (ix = 0; ix< m_imageWidth; ix++)
		for (iy = 0; iy< m_imageHeight; iy++)
		{
			if (max<m_imageRGB->GetRed(ix,iy))
				max = m_imageRGB->GetRed(ix,iy);
		}

	for (ix = 0; ix< m_imageWidth; ix++)
		for (iy = 0; iy< m_imageHeight; iy++)
		{
			s = 255.0/pow(max,p)*pow((double)m_imageRGB->GetRed(ix,iy),p);
			m_imageRGB->SetRGB(ix,iy,s,s,s);
		}
			
	Refresh();


}

void MyCanvas::Hist() //histogram equalization
{


	int num_pix = m_imageWidth*m_imageHeight;
	int ix, iy, px;
	int ss = 0; 

	temp_imageData = (unsigned char*) malloc(m_imageWidth*m_imageHeight*3);
	memcpy(temp_imageData,m_imageRGB->GetData(),m_imageWidth*m_imageHeight*3);
	temp_imageRGB = new wxImage(m_imageWidth, m_imageHeight, temp_imageData, false);	

	for (px = 0; px<=255; px++)
	{
		for (ix = 0; ix< m_imageWidth; ix++)
			for (iy = 0; iy< m_imageHeight; iy++)
				if (px == (int)temp_imageRGB->GetRed(ix,iy)) ss++;

		for (ix = 0; ix< m_imageWidth; ix++)
			for (iy = 0; iy< m_imageHeight; iy++)
			{
				if (px == (int)temp_imageRGB->GetRed(ix,iy))
					m_imageRGB->SetRGB(ix,iy,(int)(ss*255.0/num_pix+0.5),
											 (int)(ss*255.0/num_pix+0.5),
											 (int)(ss*255.0/num_pix+0.5));
			}
				
	}

	delete temp_imageRGB;
	Refresh();

}

void MyCanvas::AuThres() //automated thresholding 
{


	int num_pix_bkg, num_pix_obj, bkg_val, obj_val;
	int ix, iy;
	double thres, thres_old, epsilon;
	int ss = 0; 
	epsilon = 1;

	thres = 100;

	do {

		thres_old = thres;
		num_pix_bkg = 0;
		num_pix_obj = 0;
		bkg_val = 0;
		obj_val = 0;
		for (ix = 0; ix< m_imageWidth; ix++)
			for (iy = 0; iy< m_imageHeight; iy++)
				if (m_imageRGB->GetRed(ix,iy)<thres) {
					num_pix_obj++;
					obj_val = obj_val+ (int)(m_imageRGB->GetRed(ix,iy));
				}
				else {
					num_pix_bkg++;
					bkg_val = bkg_val+ (int)(m_imageRGB->GetRed(ix,iy));
				}

		thres = (bkg_val/(double)num_pix_bkg + obj_val/(double)num_pix_obj)/2.0;
		cout<<"thres difference: "<<abs(thres - thres_old)<<endl;		

	} while (abs(thres - thres_old)>epsilon);

	for (ix = 0; ix< m_imageWidth; ix++)
		for (iy = 0; iy< m_imageHeight; iy++)
			if (m_imageRGB->GetRed(ix,iy)<thres) {
				m_imageRGB->SetRGB(ix,iy,0,0,0);
			}
			else {
				num_pix_obj++;
				m_imageRGB->SetRGB(ix,iy,255,255,255);
			}

	Refresh();

}

void MyCanvas::Reset() //reset to inital image

{	

	memcpy(m_imageRGB->GetData(), m_imageData, m_imageWidth * m_imageHeight * 3) ;	

	Refresh(); 
}


//------------------------------------------------------------------------
void MyCanvas::BestSize() //scale window size to image size
//------------------------------------------------------------------------
{
	SetSize(m_imageWidth, m_imageHeight) ;	// ideal size for canvas
	GetParent()->SetClientSize(GetSize());	// force the main frame to show the whole canvas
}

//------------------------------------------------------------------------
void MyCanvas::OnPaint(wxPaintEvent &WXUNUSED(event))
//------------------------------------------------------------------------
// update the main window content
{

	wxPaintDC dc(this);

	if (m_imageRGB){
		m_imageBitmap = wxBitmap(*m_imageRGB, -1); 
		dc.DrawBitmap(m_imageBitmap, 0, 0) ;

	}

}

//************************************************************************
//************************************************************************
// Frame class (the main window)
//************************************************************************
//************************************************************************

//------------------------------------------------------------------------
class MyFrame: public wxFrame
//------------------------------------------------------------------------
{
public:
    MyFrame(const wxString& title, const wxPoint& pos, const wxSize& size);

protected:
	void OnQuit(wxCommandEvent& event);
    void OnAbout(wxCommandEvent& event);
	void OnOpenImage(wxCommandEvent& WXUNUSED(event) ) ;
	void OnSaveImage(wxCommandEvent & WXUNUSED(event)) ;
	void OnInvertImage(wxCommandEvent& WXUNUSED(event) ) ;
	void OnMkGrey(wxCommandEvent& WXUNUSED(event) ) ;
	void OnRescale(wxCommandEvent& WXUNUSED(event) ) ;
	void OnShift(wxCommandEvent& WXUNUSED(event) ) ;
	void OnAvergae(wxCommandEvent& WXUNUSED(event) );		
	void OnLaplacian(wxCommandEvent& WXUNUSED(event) );		
	void OnNoise_intro(wxCommandEvent& WXUNUSED(event) );	
	void OnMin(wxCommandEvent& WXUNUSED(event) );
	void OnMax(wxCommandEvent& WXUNUSED(event) );
	void OnMedian(wxCommandEvent& WXUNUSED(event) );		
	void OnNegLin(wxCommandEvent& WXUNUSED(event) );
	void OnLog(wxCommandEvent& WXUNUSED(event) );
	void OnPow(wxCommandEvent& WXUNUSED(event) );
	void OnHist(wxCommandEvent& WXUNUSED(event) );
	void OnAuThres(wxCommandEvent& WXUNUSED(event) );
	void OnReset(wxCommandEvent& WXUNUSED(event) );	
	void OnClose(wxCloseEvent& event) ;
	void OnBestSize(wxCommandEvent& WXUNUSED(event)) ;

    MyCanvas *m_canvas; // the canvas inside the main frame
	bool m_imageLoaded ;
    DECLARE_EVENT_TABLE()
};


BEGIN_EVENT_TABLE(MyFrame, wxFrame)
	EVT_MENU(ID_LOAD,  MyFrame::OnOpenImage)
	EVT_MENU(ID_SAVE,  MyFrame::OnSaveImage)
	EVT_MENU(ID_INVERT, MyFrame::OnInvertImage)
	EVT_MENU(ID_MKGREY, MyFrame::OnMkGrey)
	EVT_MENU(ID_RESCALE, MyFrame::OnRescale)
	EVT_MENU(ID_SHIFT, MyFrame::OnShift)
	EVT_MENU(ID_AVERAGE, MyFrame::OnAvergae)
	EVT_MENU(ID_LAPLACIAN, MyFrame::OnLaplacian)
	EVT_MENU(ID_NOISEINTRO, MyFrame::OnNoise_intro)
	EVT_MENU(ID_MIN,  MyFrame::OnMin)
	EVT_MENU(ID_MAX,  MyFrame::OnMin)
	EVT_MENU(ID_MEDIAN,  MyFrame::OnMedian)
	EVT_MENU(ID_NEGLIN,  MyFrame::OnNegLin)
	EVT_MENU(ID_LOG,  MyFrame::OnLog)
	EVT_MENU(ID_POW,  MyFrame::OnPow)
	EVT_MENU(ID_HIST,  MyFrame::OnHist)
	EVT_MENU(ID_AUTHRES,  MyFrame::OnAuThres)
	EVT_MENU(ID_RESET, MyFrame::OnReset)
	EVT_MENU(ID_BEST_SIZE,  MyFrame::OnBestSize)
    EVT_MENU(ID_QUIT,  MyFrame::OnQuit)
    EVT_MENU(ID_ABOUT, MyFrame::OnAbout)
	EVT_CLOSE(MyFrame::OnClose)

END_EVENT_TABLE()

//------------------------------------------------------------------------
MyFrame::MyFrame(const wxString& title, const wxPoint& pos, const wxSize& size)
       : wxFrame((wxFrame *)NULL, -1, title, pos, size)
//------------------------------------------------------------------------
{
	wxMenu *file_menu = new wxMenu();
	file_menu->Append(ID_LOAD, _T("&Open image..."));
	file_menu->Append(ID_RESET, _T("&Reset"));		
	file_menu->Append(ID_SAVE, _T("&Save image as..."));
	file_menu->Append(ID_BEST_SIZE, _T("&Best size"));
	file_menu->AppendSeparator();
	file_menu->Append(ID_ABOUT, _T("&About..."));
	file_menu->AppendSeparator();
	file_menu->Append(ID_QUIT, _T("&Exit"));

	wxMenu *file_menu2 = new wxMenu();
	file_menu2->Append(ID_MKGREY, _T("&Make Gray"));
	file_menu2->Append(ID_NOISEINTRO, _T("&Noise Introduction"));		
	file_menu2->Append(ID_INVERT, _T("&Invert image"));	
	file_menu2->Append(ID_RESCALE, _T("&Rescale"));	
	file_menu2->Append(ID_SHIFT, _T("&Shift"));	

	wxMenu *file_menu3 = new wxMenu();
	file_menu3->Append(ID_AVERAGE, _T("&Average"));	
	file_menu3->Append(ID_LAPLACIAN, _T("&Laplacian (8-neighbour)"));		
	file_menu3->Append(ID_MIN, _T("&Min Filter"));
	file_menu3->Append(ID_MAX, _T("&Max Filter"));
	file_menu3->Append(ID_MEDIAN, _T("&Median Filter"));

	wxMenu *file_menu4 = new wxMenu();
	file_menu4->Append(ID_NEGLIN, _T("&Negative Linear"));	
	file_menu4->Append(ID_LOG, _T("&Logrithm"));
	file_menu4->Append(ID_POW, _T("&Power (Gamma)"));
	file_menu4->Append(ID_HIST, _T("&Histogram Equalization"));
	file_menu4->Append(ID_AUTHRES, _T("&Automated Thresholding (Iterative Optimal)"));	

	wxMenuBar *menuBar = new wxMenuBar();

	menuBar->Append(file_menu, _T("&File"));
	menuBar->Append(file_menu2,"&Color/Dimensions");
	menuBar->Append(file_menu3,"&Filters");
	menuBar->Append(file_menu4,"&Enhancement");
	SetMenuBar( menuBar );


// create the canvas that will manage the image
	m_canvas = new MyCanvas( this, -1, wxDefaultPosition, wxDefaultSize);
	m_imageLoaded = false ;
	Centre() ;
}


//------------------------------------------------------------------------
void MyFrame::OnQuit(wxCommandEvent& WXUNUSED(event))
//------------------------------------------------------------------------
{
    Close(true) ;
}

//------------------------------------------------------------------------
void MyFrame::OnClose(wxCloseEvent& event)
//------------------------------------------------------------------------
{
	delete m_canvas ;
	event.Skip() ;
}

//------------------------------------------------------------------------
void MyFrame::OnAbout(wxCommandEvent& WXUNUSED(event))
//------------------------------------------------------------------------
{
    wxMessageBox( _T("How to \n\n- load\n- display\n- process\n- save\n\nan image with wxWidgets (2.8.7)\n\nPascal Bertolino - GIPSA-lab, Grenoble - France\npascal.bertolino@gipsa-lab.fr"),
                  _T(APP_NAME), wxOK | wxICON_INFORMATION ) ;
}


//------------------------------------------------------------------------
void MyFrame::OnInvertImage(wxCommandEvent& WXUNUSED(event))
//------------------------------------------------------------------------
{
	if (m_imageLoaded)
	    m_canvas->InvertImage() ;
}

void MyFrame::OnMkGrey(wxCommandEvent& WXUNUSED(event))
//------------------------------------------------------------------------
{
	if (m_imageLoaded)
	    m_canvas->MkGrey() ;
}

//------------------------------------------------------------------------
void MyFrame::OnRescale(wxCommandEvent& WXUNUSED(event))
//------------------------------------------------------------------------
{
	if (m_imageLoaded)
	    m_canvas->Rescale() ;
}

void MyFrame::OnShift(wxCommandEvent& WXUNUSED(event))
//------------------------------------------------------------------------
{
	if (m_imageLoaded){


		wxString shift_am = wxGetTextFromUser("Please enter the shift mmount:","Image Shift","0");

		wxString number(shift_am);
		double value;
		if(number.ToDouble(&value))
			m_canvas->Shift(value) ;


	}
}

void MyFrame::OnAvergae(wxCommandEvent& WXUNUSED(event))
{
	if (m_imageLoaded)
	    m_canvas->Average() ;	
}

void MyFrame::OnLaplacian(wxCommandEvent& WXUNUSED(event))
{
	if (m_imageLoaded)
	    m_canvas->Laplacian() ;	
}

void MyFrame::OnNoise_intro(wxCommandEvent& WXUNUSED(event))
//------------------------------------------------------------------------
{
	if (m_imageLoaded){


		wxString noise_am = wxGetTextFromUser("Please enter the noise amount (number of random iteration):","num of iteration","0");

		wxString number(noise_am);
		long num;
		if(number.ToLong(&num))
			m_canvas->Noise_intro(num) ;


	}
}

void MyFrame::OnMin(wxCommandEvent& WXUNUSED(event))
{
	if (m_imageLoaded)
	    m_canvas->Min() ;	
}

void MyFrame::OnMax(wxCommandEvent& WXUNUSED(event))
{
	if (m_imageLoaded)
	    m_canvas->Max() ;	
}

void MyFrame::OnMedian(wxCommandEvent& WXUNUSED(event))
{
	if (m_imageLoaded)
	    m_canvas->Median() ;	
}

void MyFrame::OnNegLin(wxCommandEvent& WXUNUSED(event))
{
	if (m_imageLoaded)
	    m_canvas->NegLin() ;	
}

void MyFrame::OnLog(wxCommandEvent& WXUNUSED(event))
{
	if (m_imageLoaded)
	    m_canvas->Log() ;	
}

void MyFrame::OnPow(wxCommandEvent& WXUNUSED(event))
{
	if (m_imageLoaded)
	    m_canvas->Pow() ;	
}

void MyFrame::OnHist(wxCommandEvent& WXUNUSED(event))
{
	if (m_imageLoaded)
	    m_canvas->Hist() ;	
}

void MyFrame::OnAuThres(wxCommandEvent& WXUNUSED(event))
{
	if (m_imageLoaded)
	    m_canvas->AuThres() ;	
}

void MyFrame::OnReset(wxCommandEvent& WXUNUSED(event))
{
	if (m_imageLoaded)
	    m_canvas->Reset() ;	
}


//------------------------------------------------------------------------
void MyFrame::OnOpenImage(wxCommandEvent& WXUNUSED(event) )
//------------------------------------------------------------------------
{
	wxBitmap bitmap;

	wxString filename = wxFileSelector("Select file",_T(""),_T(""),_T(""), _T("All files (*.*)|*.*") );
	if ( !filename.empty() )
	{
		m_canvas->LoadImage(filename) ;
		m_imageLoaded = true ;
	}
}

//------------------------------------------------------------------------
void MyFrame::OnSaveImage(wxCommandEvent & WXUNUSED(event))
//------------------------------------------------------------------------
{
//	char str[128] = "" ; // proposed file name

	if (!m_imageLoaded)
		return ;

	wxString filename = wxFileSelector(_T("Save image as"),_T(""),_T(""),_T("*.bmp"), _T("BMP files (*.bmp)|*.bmp|GIF files (*gif)|*.gif|JPEG files (*jpg)|*.jpg|PNG files (*png)|*.png|TIFF files (*tif)|*.tif|XPM files (*xpm)|*.xpm|All files (*.*)|*.*"), wxFD_SAVE);
	if ( !filename.empty() )
		m_canvas->SaveImage(filename) ;
}

//------------------------------------------------------------------------
void MyFrame::OnBestSize(wxCommandEvent& WXUNUSED(event))
//------------------------------------------------------------------------
{
    m_canvas->BestSize() ;
}

//************************************************************************
//************************************************************************
// Application class
//************************************************************************
//************************************************************************

//------------------------------------------------------------------------
class MyApp: public wxApp
//------------------------------------------------------------------------
{
    virtual bool OnInit() ;
};

IMPLEMENT_APP(MyApp) // macro that contains the main() function


//------------------------------------------------------------------------
bool MyApp::OnInit()
//------------------------------------------------------------------------
{
//support all available image formats
	wxInitAllImageHandlers() ;

    MyFrame *frame = new MyFrame(_T(APP_NAME), wxDefaultPosition, wxSize(400,300)) ;
    frame->CreateStatusBar();
    frame->SetStatusText("'Education is what remains after one has forgotten what one has learned in school.'--A.Einstein");

    frame->Show(true) ;
    SetTopWindow(frame) ;

    return true ;
}
