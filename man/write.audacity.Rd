\name{write.audacity}

\alias{write.audacity}

\title{Audacity audio markers export}

\description{Write audio markers to be imported by Audacity.}

\usage{write.audacity(x, filename)}

\arguments{
  \item{x}{a data frame with the three or five colums, see details.}
  \item{filename}{name of the \code{.txt} file. (by default the name of
    \code{x}).}
}

\details{The input \code{x} object should be a data frame with two or
  three columns depending on whether the markers include frequency
  limits or not : \cr
  \itemize{
    \item{time limits only:
      \enumerate{
	\item{text label of each marker,}
	\item{time marker of the beginning of each marker,}
	\item{time marker of the end of each marker.}
      }
    }
    \item{time and frequency limits:
      \enumerate{
	\item{text label of each marker,}
	\item{time marker of the beginning of each marker,}
	\item{time marker of the end of each marker,}
        \item{lower frequency limit of each marker,}
	\item{higher frequency limit of each marker.}
      }
    }
  }
}


\value{A \code{.txt} file is generated to be imported as a markers in Audacity.}

\note{Naming the columns of \code{x} is not necessary.}

\references{Audacity is a free software
  distributed under the terms of the GNU General Public License.\cr
  Web site: \url{https://audacityteam.org/. }
}

\author{Jerome Sueur}

\seealso{\code{\link{read.audacity}}}

\examples{
## 3 markers, time only
t1 <- c(9.2, 16.2, 24.4)
t2 <- c(11.7, 18.7, 26.9)
label <- c("a", "b", "c")
df <- data.frame(label, t1, t2)
write.audacity(df, filename="test-time.txt")
## 3 markers, time and frequency
t1 <- c(9.4, 15.2, 24.9)
t2 <- c(10.54, 16.6, 26.1)
f1 <- c(1703.4, 3406.8, 1608.8)
f2 <- c(7476.2, 8517.2, 5110.3)
label <- c("a", "b", "c")
dff <- data.frame(label, t1, t2, f1, f2)
write.audacity(dff, filename="test-time-frequency.txt")
## delete files
unlink(c("test-time.txt", "test-time-frequency.txt"))
}

\keyword{input}
\keyword{data}
