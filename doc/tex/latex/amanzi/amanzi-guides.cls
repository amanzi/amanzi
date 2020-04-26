% ==========================================================================
%
%  Amanzi: Class for the Theory Guide and Developers Guide 
%
% ==========================================================================

\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{guides}[2020/03/28 Guides class for ASCEM]

% based on srs.cls for docs repository

% ==========================================================================
% -------------------------------------------------
%  Options that we don't inherit from article.cls
%
%  - note this means any functionality associated 
%    with these options in article.cls will need 
%    to be provided by this class
% -------------------------------------------------

\newif{\if@draft}
\@draftfalse
\DeclareOption{draft}{\@drafttrue}

\newif{\if@final}
\@finalfalse
\DeclareOption{final}{\@finaltrue}

\newif{\if@sidenotes}
\@sidenotesfalse
\DeclareOption{sidenotes}{\@sidenotestrue}

%
% [final] from article.cls defines
%
\setlength\overfullrule{0pt}

% -------------------------------------------------
%  Pass all other options through to article
% -------------------------------------------------

\DeclareOption*{\PassOptionsToClass{\CurrentOption}{article}}

% -------------------------------------------------
%  Process these options 
% -------------------------------------------------

\ProcessOptions\relax

% -------------------------------------------------
%  Load the base class
% -------------------------------------------------

\if@sidenotes
  \LoadClass[11pt]{article}
\else
  \LoadClass[12pt]{article}
\fi

% ==========================================================================
% -------------------------------------------------
%  Packages
% -------------------------------------------------

% ---------------------
%  Handy Stuff
% ---------------------

\RequirePackage{ifthen}
\RequirePackage{ifpdf}
\RequirePackage{verbatim}
\RequirePackage{xspace}

% ---------------------
%  Fonts
% ---------------------

\RequirePackage{times}
\RequirePackage{amsfonts,amsmath,amssymb}

% ---------------------
%  Geometry
% ---------------------

\RequirePackage{geometry}

\if@sidenotes
\geometry{pdftex, includemp=false, portrait, 
          headsep=24pt, headheight=15pt, footskip=27pt,
          marginparwidth=2.5in,
          top=1.00in, bottom=1.00in,
          left=1.00in, right=3.00in}
\else
\geometry{pdftex, includemp=false, portrait, 
          headsep=24pt, headheight=15pt, footskip=27pt,
          top=1.00in, bottom=1.00in,
          left=1.00in, right=1.00in}
\fi

% ---------------------
%  Color and Graphics
% ---------------------

\ifpdf
  \PassOptionsToPackage{pdftex}{xcolor}
  \PassOptionsToPackage{pdftex,final}{graphicx}
\else
  \PassOptionsToPackage{dvips}{xcolor}
  \PassOptionsToPackage{dvips,final}{graphicx}
\fi

\RequirePackage{xcolor}
\RequirePackage{graphicx}

% ---------------------
%  Hyperlinks
% ---------------------

% Note that the colornames we use here are provided 
% by the options we use in the xcolor package (above).
% - colornames failed for some people!

% MidnightBlue
\definecolor{ASCEMb}{rgb}{.10,.10,.50}
% DarkOrchid/Dark Magenta
\definecolor{ASCEMp}{rgb}{.55, 0, .55}
% SeaGreen
\definecolor{ASCEMg}{rgb}{.18,.55,.34}

\PassOptionsToPackage{
  colorlinks=true,
  linkcolor=ASCEMb,
  urlcolor=ASCEMp,
  citecolor=ASCEMg,
}{hyperref}

\ifpdf
  \PassOptionsToPackage{pdftex,draft=false}{hyperref}
\else
  \PassOptionsToPackage{dvips,draft=false}{hyperref}
\fi
\RequirePackage{hyperref}

% ---------------------
%  Include pages 
% ---------------------

\RequirePackage{pdfpages}     

% ---------------------
%  ToDo notes
% ---------------------

\if@sidenotes
  \PassOptionsToPackage{textsize=small,colorinlistoftodos}{todonotes}
\else
  \PassOptionsToPackage{textsize=small,disable}{todonotes}
\fi

\RequirePackage{todonotes}

% ---------------------
%  Draft WaterMark
% ---------------------

\if@draft
  \RequirePackage{draftwatermark}
  \SetWatermarkLightness{0.95}
  \SetWatermarkScale{1.1}
  %
  \SetWatermarkText{%
    \parbox[t]{7in}{%
      DRAFT\\[-2in]
      \hspace*{0.25in}
      {\Huge Amanzi Team Internal Distribution Only}
    }
  }
\fi

% ---------------------
%  Equation numbers:
% ---------------------

\numberwithin{equation}{section}

% ---------------------
%  Return Booleans:
% ---------------------

\newboolean{sidenotes}
\if@sidenotes
  \setboolean{sidenotes}{true}
\else
  \setboolean{sidenotes}{false}
\fi

% ==========================================================================
