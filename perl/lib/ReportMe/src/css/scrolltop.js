$(document).ready(function () {
  $("#goTop").hide();

  $(window).scroll(function () {
    if ($(window).scrollTop() > 0) {
      $("#goTop").fadeIn(800);
    } else {
      $("#goTop").fadeOut(800);
    }
  });

  $("#goTop").click(function () {
    $("body,html").animate({ scrollTop: 0 }, 500);
  });
});
