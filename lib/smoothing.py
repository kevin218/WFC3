import scipy.signal as sps
import numpy as np

def smoothing(im, nPix, sPix, mask=np.array([]), gk=np.array([]), mirror=False):
	"""
	This function performs 2D smoothing by convolving a masked image with a kernel.

	Parameters
	----------
	im :	2D array
			Image to be smoothed
	nx,ny : Int
			Width of smoothing function along x & y, in pixels
	sx,sy : Int
			Std dev of Gaussian kernel along x & y, in pixels
	mask :  2D array
			1s are used in smoothing
	gk :	4D array, (im.shape, nx, ny)
			Guassian-distributed kernel for each image position
	mirror: boolean
			Set to True to have the padded image values mirror the existing values
	
	Returns
	-------
	output :	2D array
				Smoothed image

	Examples
	--------
	None

	Revisions
	---------
	2010-06-19  Kevin Stevenson, UCF
				kevin218@knights.ucf.edu
				Original version
	"""
	
	ny,nx = nPix
	
	if sPix is not None:
		[sy,sx] = sPix
	else:
		[sy,sx] = 1,1
	
	sizey, sizex = im.shape
	if (mask.size == 0) and (gk.size == 0):
		#Perform normal 2D concolution
		gk = gauss_kernel((ny,nx), (sy,sx))
		return sps.convolve2d(im, gk, mode='same')
	'''
	#FINDME: function runs out of memory with large arrays
	elif gk.size == 0:
		#Calculate kernel since none given
		gk = gauss_kernel_mask((ny,nx), (sy,sx), mask)
	'''
	#Pad image initially with zeros
	newim = np.zeros(np.array(im.shape) + 2*np.array((ny, nx)))
	newim[ny:-ny, nx:-nx] = im
	if mirror == True:
		#Pad image with mirrored data
		newim[   :  ,   :nx] = np.fliplr(newim[	 :	,   nx:2*nx])
		newim[   :  ,-nx:  ] = np.fliplr(newim[	 :	,-2*nx: -nx])
		newim[   :ny,   :  ] = np.flipud(newim[   ny:2*ny,	 :	])
		newim[-ny:  ,   :  ] = np.flipud(newim[-2*ny: -ny,	 :	])
	#Perform 2D convolution where mask = 1
	smooth = np.zeros((sizey,sizex))
	for j in range(sizey):
		print(j,sizey)
		for i in range(sizex):
			if mask[j,i] == 1:
				#Calculate kernel since none given
				gk = gauss_kernel_mask2((ny,nx), (sy,sx), (j,i), mask)
				smooth[j,i] = np.sum(gk * newim[j:j+2*ny+1, i:i+2*nx+1])
	
	return smooth

def gauss_kernel(nPix, sPix):
	'''
		nPix = [ny, nx] # number of pixels in y and x direction
		sPix = [sy, sx] # widths in pixels in y and x direction
	'''
	
	ny,nx = nPix
	sy,sx = sPix
	
	x, y  = np.mgrid[-nx:nx+1, -ny:ny+1]
	gauss = np.exp(-(0.5*(x/sx)**2 + 0.5*(y/sy)**2))
	return gauss/gauss.sum()

def gauss_kernel_mask2(nPix, sPix, j_i_inds, mask):
	'''
		nPix 	 = [ny, nx] # number of pixels in y and x direction
		sPix 	 = [sy, sx] # widths in pixels in y and x direction
		j_i_inds = [j,i]    # indices to iterate over with the "newmask"
	'''
	ny,nx   = nPix
	sy,sx   = sPix
	j,i		= j_i_inds
	
	#x,	 y	 = 1.*np.mgrid[-nx:nx+1, -ny:ny+1]
	y,	 x	 = 1.*np.mgrid[-ny:ny+1, -nx:nx+1]
	#Create new mask with padding
	sizey, sizex = mask.shape
	newmask	  = np.zeros(np.array((sizey, sizex)) + 2*np.array((ny, nx)))
	#Copy the mask into the middle of the new mask, leaving padding as zeros
	newmask[ny:-ny, nx:-nx] = mask
	gauss = np.exp(-(0.5*(x/sx)**2 + 0.5*(y/sy)**2)) * newmask[j:j+2*ny+1, i:i+2*nx+1]
	gsum  = gauss.sum()
	if gsum > 0:
		kernel = gauss/gsum
	else:
		kernel = 0.
		#kernel = gauss
	
	return kernel

def gauss_kernel_mask(nPix, sPix, mask):
	'''
		nPix = [ny, nx] # number of pixels in y and x direction
		sPix = [sy, sx] # widths in pixels in y and x direction
	'''
	
	ny,nx = nPix
	sy,sx = sPix
	
	sizey, sizex = mask.shape
	x,	 y	 = 1.*np.mgrid[-nx:nx+1, -ny:ny+1]
	kernel	   = np.zeros((sizey, sizex, 2*ny+1, 2*nx+1))
	if np.size(sy) == 1:
		sy = sy*np.ones((sizey, sizex))
		sx = sx*np.ones((sizey, sizex))
	#Create new mask with padding
	newmask	  = np.zeros(np.array((sizey, sizex)) + 2*np.array((ny, nx)))
	#Copy the mask into the middle of the new mask, leaving padding as zeros
	newmask[ny:-ny, nx:-nx] = mask
	for j in range(sizey):
		for i in range(sizex):
			gauss = np.exp(-(0.5*(x/sx[j,i])**2 + 0.5*(y/sy[j,i])**2)) * newmask[j:j+2*ny+1, i:i+2*nx+1]
			gsum  = gauss.sum()
			if gsum > 0:
				kernel[j,i] = gauss/gsum
			else:
				kernel[j,i] = gauss
	
	return kernel

