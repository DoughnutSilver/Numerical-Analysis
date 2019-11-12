using VideoIO
using Makie
#viewcam()      #顔が映る
#VideoIO.CAMERA_DEVICESD #camera avalable

f=opencamera()
img=read(f)     #pass RGB array
close(f)
